# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:45:52 2024

@author: mre283
"""

import pandas as pd
import cbmpy
import numpy as np
import scipy as sp
import sympy
import matplotlib.pyplot as plt

class separate_cat_ana():
    def __init__(self, template_file_path, num_check=True, all_cofactors=True, multiple_EFMs=False):
        """
        Initialize the class with a given template file path and optional parameters.
        
        :param template_file_path: Path to the Excel file containing the template data
        :param num_check: Boolean to indicate if numerical checks should be performed
        :param all_cofactors: Boolean to indicate if all cofactors should be included
        :param multiple_EFMs: Boolean to indicate if multiple Elementary Flux Modes (EFMs) are to be considered
        """
        
        # Read the Excel file without header and setting the first column as index
        readfile = pd.read_excel(template_file_path, header=None, index_col=0)
        self.num_check = num_check
        self.all_cofactors = all_cofactors
        
        # Extracting various parameters from the Excel file
        self.simtitle = readfile.loc['simulation title', 1]
        self.modelname = readfile.loc['model name', 1]
        self.modelpath = readfile.loc['model file location', 1]
        self.reversiblereacs = readfile.loc['reversible reactions present', 1]
        self.objective_flux = readfile.loc['objective flux', 1]
        self.optim_sense = readfile.loc['optimization sense', 1]
        self.biomass_reac = readfile.loc['biomass reaction', 1]
        self.multiple_EFMs = multiple_EFMs
        
        # Extracting and processing constraints
        constraints = readfile.loc['constraints', 1].split('; ')
        self.constraints_dict = {}
        for i in constraints:
            temp = i.split(', ')
            self.constraints_dict[temp[0]] = (float(temp[1]), float(temp[2]))
        
        # Extracting and processing energy currencies
        energy_currencies = readfile.loc['Energy currencies', 1].split('\n')
        self.currencies_dict = {}
        for curr in energy_currencies:
            lcurr = curr.split(':')
            tuples = lcurr[1][1:].split('), (')
            
            tuplist = []
            for tp in tuples:
                tp = tp.replace('[(', '')
                tp = tp.replace(')]', '')
                els = tp.split(', ')
                tuplist.append((els[0], float(els[1])))
            
            self.currencies_dict[lcurr[0]] = tuplist
        
        # Extracting macromolecule-related data
        self.split_macromolecules = readfile.loc['Split macromolecules', 1]
        if self.split_macromolecules is True:
            self.macromolecule_exch_ids = readfile.loc['Macromolecule exchange reactions', 1].split(', ')
            self.macromolecule_species_ids = readfile.loc['Macromolecule species ids', 1].split(', ')
            GAM = readfile.loc['Growth associated maintenance reaction', 1].split(', ')
            self.GAM_id, self.GAM_value = GAM[0], np.float64(GAM[1])
        
        # Step one: load the model and set constraints
        self.model = cbmpy.loadModel(self.modelpath)  # Load the metabolic model
    
        # Set reversible reactions if specified
        if self.reversiblereacs is True:
            if self.objective_flux != self.biomass_reac:
                self.model.getReaction(self.objective_flux).reversible = True
    
        # Adjust the objective flux based on the optimization sense
        if self.model.getReaction(self.objective_flux).reversible is True:
            if self.optim_sense == 'minimize':
                self.objective_flux = self.objective_flux + '_fwd'
            elif self.optim_sense == 'maximize':
                self.objective_flux = self.objective_flux + '_rev'
                self.optim_sense = 'minimize'
            else:
                raise Exception('Optimization sense should be either maximize or minimize')
    
        # Initialize empty lists for ATP and cofactor dataframes and required reactions
        self.atp_df = []
        self.cofactor_df = []
        self.required = []
        
        # Set reaction bounds based on constraints
        for const in list(self.constraints_dict.keys()):
            self.model.setReactionBounds(const, self.constraints_dict[const][0], self.constraints_dict[const][1])


    def do_cat_ana_split(self):
        """
        Performs splitting of catabolism and anabolism according to the options supplied to this class.
        
        Steps:
        1. Analyze the model to get active fluxes and constraints.
        2. Split the EFMs if multiple EFMs are required.
        3. For each model, split catabolism and anabolism for all cofactors using the perform_splitting_multiplecofactors function.
        4. If only required cofactors are needed, reduce the model accordingly.
        5. If desired, split the model into macromolecule production and split each macromolecule into catabolism and anabolism.
    
        Returns:
        None
        """
    
        # Step 1: Analyze the model to get active fluxes, the number of active fluxes, and active constraints
        active_fluxdict, N_active, active_constraints = self.analyze_model()
        
        # Initialize dictionaries to store results for different models and dataframes
        model_dict = {}
        df_dict = {}
        df_atp_dict = {}
        df_cofactor_dict = {}
    
        # Step 2: Handle multiple EFMs if specified
        if len(list(active_constraints.keys())) > 1:
            if not self.multiple_EFMs:
                print(active_constraints)
                raise Exception('More than 1 EFM found')
            
            unbounded_dict = {}
            df_EFMfluxes = pd.DataFrame(np.zeros(len(active_fluxdict.keys())), columns=[list(active_constraints.keys())[0]], index=list(active_fluxdict.keys()))
    
            # Process each active constraint to generate EFMs
            for k, v in active_constraints.items():
                unboundeds = [l for l in list(active_constraints.keys()) if l != k]
                unbounded_dict[k] = unboundeds
                EFM = self.generate_EFM(self.active_model, k, [v, v], unboundeds, self.objective_flux, [0, 1000], self.optim_sense)
                EFMfluxes = EFM.getReactionValues()
                df_EFMfluxes.loc[:, k] = EFMfluxes
                EFM = self.get_activemodel(EFM)
                model_dict[k] = EFM
                for r in unboundeds:
                    EFM.setReactionBounds(r, 0, 1000)
            
            # Identify zero flux EFMs
            df_EFMfluxes_zero = df_EFMfluxes[(np.abs(df_EFMfluxes) < 1e-7).any(axis=1)]
            
            # Scale the fluxes
            scale_fluxes = {}
            for EFM in df_EFMfluxes_zero.columns:
                for v in df_EFMfluxes_zero.index:
                    if np.abs(df_EFMfluxes_zero.loc[v, EFM]) > 1e-7 and df_EFMfluxes_zero.loc[v].astype(bool).sum() == 1:
                        scale_fluxes[EFM] = v
            
            df_EFMfluxes_scaled = df_EFMfluxes
            array_EFMs = np.array(df_EFMfluxes)
            array_orgflux = np.array(list(active_fluxdict.values()))
            fractions = np.linalg.pinv(array_EFMs) @ array_orgflux
            
            for i, f in enumerate(fractions):
                df_EFMfluxes_scaled.iloc[:, i] = f * df_EFMfluxes.iloc[:, i]
        
        # If there is only 1 EFM active, the active model is taken
        else:
            model_dict['tot'] = self.active_model.clone()
        
        # Step 3: Perform splitting for each model (each model corresponds to an EFM)
        biomass_analysis_obj_dict = {}
        for name, m in model_dict.items():
            if len(list(model_dict.keys())) > 1:
                active_fluxdict = df_EFMfluxes_scaled[name].to_dict()
            
            df_biomass, df_ATP, df_cofactors, biomass_analysis_obj = self.perform_splitting_multiplecofactors(m, active_fluxdict, self.currencies_dict)
            df_dict[name] = df_biomass
            df_atp_dict[name] = df_ATP
            df_cofactor_dict[name] = df_cofactors
            biomass_analysis_obj_dict[name] = biomass_analysis_obj
    
        self.atp_df = []
        self.cofactor_df = []
        self.result_df = []
        
        if self.all_cofactors:
            if len(list(model_dict.keys())) > 1:
                # Step 4: Combine results for multiple EFMs
                df_tot = pd.DataFrame(list(active_fluxdict.values()), index=list(active_fluxdict.keys()), columns=['org fluxdict'])
                org_cols = list(df_dict.values())[0].columns
    
                # Concatenate results from each EFM
                for name, df_EFM in df_dict.items():
                    df_EFM = df_EFM.rename(columns={i: f'{i}_{name}' for i in df_EFM.columns})
                    df_tot = pd.concat([df_tot, df_EFM], axis=1).fillna(0)
    
                # Sum the results for each original column
                for col in org_cols:
                    df_tot[f'{col}_tot'] = sum(df_tot[f'{col}_{name}'] for name in df_dict)
                df_tot['biomass sumcheck'] = df_tot['biomass_tot'] - df_tot['org fluxdict']
                df_dict['tot'] = df_tot
                self.result_df = df_tot
    
                # Concatenate and sum results for cofactors and ATP
                df_cofactor_tot = self.concat_and_sum(df_cofactor_dict, org_cols)
                df_cofactor_dict['tot'] = df_cofactor_tot
                self.cofactor_df = df_cofactor_tot
                
                df_atp_tot = self.concat_and_sum(df_atp_dict, org_cols)
                df_atp_dict['tot'] = df_atp_tot
                self.atp_df = df_atp_tot
            
            else:
                self.result_df = list(df_dict.values())[0]
                self.atp_df = list(df_atp_dict.values())[0]
                self.cofactor_df = list(df_cofactor_dict.values())[0]
        
        else:
            # Step 5: Perform Flux Variability Analysis (FVA) to identify required cofactors
            if len(list(model_dict.keys())) > 1:
                df_tot, df_tot_atp, df_tot_cofactor = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
                result_org_cols, atp_org_cols, cofactor_org_cols = set(), set(), set()
                df_dict2 = {}
    
                # Define required exchanges for each EFM
                for name, model in model_dict.items():
                    result_df, atp_df, cofactor_df, required = self.define_required_exchanges(
                        df_dict[name]['biomass'].to_dict(), 
                        biomass_analysis_obj_dict[name], 
                        df_dict[name], 
                        self.currencies_dict
                    )
                    df_dict2[name] = [result_df, atp_df, cofactor_df]
                    result_org_cols.update(result_df.columns)
                    atp_org_cols.update(atp_df.columns)
                    cofactor_org_cols.update(cofactor_df.columns)
                    
                # Concatenate results for each EFM
                for name in model_dict:
                    result_df, atp_df, cofactor_df = df_dict2[name]
                    df_tot = self.concat_results(df_tot, result_df, name)
                    df_tot_atp = self.concat_results(df_tot_atp, atp_df, name)
                    df_tot_cofactor = self.concat_results(df_tot_cofactor, cofactor_df, name)
    
                # Sum the results for each original column
                self.result_df = self.sum_columns(df_tot, result_org_cols, model_dict)
                self.atp_df = self.sum_columns(df_tot_atp, atp_org_cols, model_dict)
                self.cofactor_df = self.sum_columns(df_tot_cofactor, cofactor_org_cols, model_dict)
            
            else:
                self.result_df, self.atp_df, self.cofactor_df, self.required = self.define_required_exchanges(
                    active_fluxdict, 
                    biomass_analysis_obj, 
                    list(df_dict.values())[0], 
                    self.currencies_dict
                )
        
        self.macromolecule_dfs = []
    
        # Step 6: Split macromolecules if required
        if self.split_macromolecules:
            if len(list(model_dict.keys())) > 1:
                raise Exception('Macromolecule splitting with more than one EFM is not supported')
            
            if self.all_cofactors:
                self.macromolecule_dfs, self.atp_macromolecule_dfs, self.cofactor_macromolecule_dfs = self.do_split_macromolecules(self.result_df, self.currencies_dict)
            else:
                self.macromolecule_dfs, self.atp_macromolecule_dfs, self.cofactor_macromolecule_dfs = self.do_split_macromolecules(self.result_df, self.required)

    def concat_and_sum(self, df_dict, org_cols):
        """
        Helper function to concatenate and sum columns for multiple EFMs.
        """
        df_tot = pd.DataFrame()
        for name, df_EFM in df_dict.items():
            df_EFM = df_EFM.rename(columns={i: f'{i}_{name}' for i in df_EFM.columns})
            df_tot = pd.concat([df_tot, df_EFM], axis=1).fillna(0)
        
        for col in org_cols:
            df_tot[f'{col}_tot'] = sum(df_tot[f'{col}_{name}'] for name in df_dict)
        
        return df_tot
    
    def concat_results(self, df_tot, result_df, name):
        """
        Helper function to concatenate results.
        """
        result_df = result_df.rename(columns={i: f'{i}_{name}' for i in result_df.columns})
        return pd.concat([df_tot, result_df], axis=1).fillna(0)
    
    def sum_columns(self, df_tot, org_cols, model_dict):
        """
        Helper function to sum columns across multiple EFMs.
        """
        for col in org_cols:
            df_tot[f'{col}_tot'] = sum(df_tot[f'{col}_{name}'] for name in model_dict)
        return df_tot
    
    def generate_dfs(self):
        """
        Returns the dataframes created in do_cat_ana_split in a list. 
        If the macromolecules are not split, the list is: 
            [complete flux distribution, summed ATP fluxes, summed NAD(P)H fluxes]
        If the macromolecules are split, the list is: 
            [complete flux distribution, summed ATP fluxes, summed NAD(P)H fluxes, 
             list of macromolecule flux distributions, list of summed ATP fluxes for each macromolecule, 
             list of summed NAD(P)H fluxes for each macromolecule]
        The lists of the macromolecules are in the order: 
            protein, lipid, lps, rna, dna, murein, ion, cofactor, GAM
    
        Returns
        -------
        list
            A list of pandas.DataFrames as described above.
        """
        
        # If macromolecules are not split, return the main result dataframes
        if not self.split_macromolecules:
            return [self.result_df, self.atp_df, self.cofactor_df]
        else:
            # If macromolecules are split, return the main result dataframes along with the macromolecule-specific dataframes
            return [self.result_df, self.atp_df, self.cofactor_df, 
                    self.macromolecule_dfs, self.atp_macromolecule_dfs, self.cofactor_macromolecule_dfs]
        
        
        
    def generate_MCEQs(self, df, columns=None, biomass_reac='R_EX_BIOMASS', additional_active_metabolites=[]):
        """
        Generates macrochemical equations (MCEQs) for a dataframe of flux distributions.
    
        This method computes the macrochemical equations for each flux distribution provided in the dataframe.
        The MCEQs are derived using the flux distributions, a specified biomass reaction, and any additional 
        active metabolites provided. Each MCEQ represents a macrochemical balance for a specific flux distribution.
    
        Parameters
        ----------
        df : pandas.DataFrame
            A dataframe where each column represents a different flux distribution.
            
        columns : list of str, optional
            A list of column names to include from the dataframe. If None, all columns are used. The default is None.
            
        biomass_reac : str, optional
            The identifier for the biomass reaction used in generating MCEQs. The default is 'R_EX_BIOMASS'.
            
        additional_active_metabolites : list of str, optional
            A list of additional active metabolites to consider in the MCEQs. The default is an empty list.
    
        Returns
        -------
        MCEQ_dict : dict
            A dictionary where keys are column names from the dataframe and values are the corresponding macrochemical equations.
        """
        
        # Use all columns if none are specified
        if columns is None:
            columns = df.columns
        
        MCEQ_dict = {}
        
        for c in columns:
            # Calculate the macrochemical equation for the current flux distribution
            mceq = self.get_MCEQ(
                self.active_model, 
                biomass_reac=biomass_reac, 
                additional_active_metabolites=additional_active_metabolites, 
                flux_dist=df[c]
            )
            # Store the result in the dictionary with the column name as the key
            MCEQ_dict[c] = mceq
            
        return MCEQ_dict

    def analyze_model(self):
        """
        Analyzes the metabolic model by performing various steps to optimize and assess its performance.
    
        This method performs the following steps:
        1. Splits reversible reactions into two irreversible reactions if necessary.
        2. Sets the objective flux and optimization sense.
        3. Runs the model using the simplex method to obtain a parsimonious solution.
        4. Retrieves active fluxes and constraints.
        5. Removes inactive fluxes to create an active model.
    
        Returns
        -------
        active_fluxdict : dict
            A dictionary of active fluxes where the keys are the flux identifiers and the values are the flux rates.
        N_active : int
            The number of active fluxes in the model.
        active_constraints : list
            A list of active constraints in the model.
        """
        
        # Step 1b: Split all reversible reactions into two irreversible ones if the model contains reversible reactions
        if self.reversiblereacs:
            self.model = cbmpy.CBTools.splitReversibleReactions(self.model)
        
        # Set the objective flux and the optimization sense (e.g., maximization or minimization)
        self.model.setObjectiveFlux(self.objective_flux, osense=self.optim_sense)
        
        # Step 2: Run the model using the simplex method to obtain a parsimonious solution
        result = cbmpy.CBGLPK.glpk_analyzeModel(self.model, method='s')
        
        # Retrieve the number of active fluxes in the model
        N_active = self.get_activeN(self.model)
        
        # Retrieve a dictionary of active fluxes and their rates
        active_fluxdict = self.get_activeFluxdict(self.model)
        
        # Perform a rank check on the active fluxes and constraints
        self.get_rankcheck(self.model, active_fluxdict, N_active, self.num_check)
        
        # Retrieve a list of active constraints in the model
        active_constraints = self.get_activeconstraints(self.model, active_fluxdict)
        
        # Step 3: Remove all inactive fluxes to create an active model
        self.active_model = self.get_activemodel(self.model)
        
        return active_fluxdict, N_active, active_constraints


    def perform_splitting_multiplecofactors(self, model, active_fluxdict, currencies):
        """
        Performs the splitting of reactions involving multiple cofactors in the given model.
    
        This method performs a series of steps to analyze and modify the metabolic model with respect to
        the biomass reaction and cofactor currencies. It includes splitting reactions based on cofactors
        and calculating associated metrics.
    
        Parameters
        ----------
        model : object
            The metabolic model to be analyzed and modified.
            
        active_fluxdict : dict
            A dictionary of active fluxes where the keys are flux identifiers and the values are flux rates.
            
        currencies : dict
            A dictionary of cofactor currencies where keys are cofactor identifiers and values are associated data.
    
        Returns
        -------
        df_biomass : pandas.DataFrame
            A dataframe with biomass flux values, indexed by flux identifiers.
            
        df_ATP : pandas.DataFrame
            A dataframe with calculated ATP-related data.
            
        df_cofactors : pandas.DataFrame
            A dataframe with calculated cofactor-related data.
            
        biomass_analysis : macromolecule_analysis
            The analysis object used to perform the macromolecule-related computations and modifications.
        """
        
        # Create an analysis object for the biomass reaction
        biomass_analysis = macromolecule_analysis(model, 'biomass', biomass_reac=self.biomass_reac)
        
        # Print the current biomass flux value
        print({self.biomass_reac: active_fluxdict[self.biomass_reac]})
        
        # Perform analysis on the biomass reaction with the provided flux value
        biomass_analysis.perform_analysis(
            {self.biomass_reac: active_fluxdict[self.biomass_reac]}, 
            optimize_flux=self.objective_flux
        )
        
        # Create a dataframe for biomass flux values
        df_biomass = pd.DataFrame(columns=['biomass'], index=active_fluxdict.keys())
        df_biomass['biomass'] = active_fluxdict.values()
        biomass_analysis.df = df_biomass
        
        # Build the stoichiometric matrix from the active model and get the row indices
        metabs = biomass_analysis.active_model.buildStoichMatrix(matrix_type='numpy', only_return=True).row
        
        # Filter out currencies that are not in the metabolite list
        currencies_temp = {}
        for k, v in currencies.items():
            for m in v:
                if m[0] not in metabs:
                    currencies_temp[k] = v
                    break
        
        # Remove filtered currencies from the original dictionary
        for k in list(currencies_temp.keys()):
            currencies.pop(k)
        
        count = 0
        booll = False
        
        # Iterate over remaining currencies to perform analysis
        for k, v in currencies.items():
            count += 1
            if count >= len(currencies):
                booll = True
            
            # Print the current model
            print(model)
            
            # Perform separation and analysis based on cofactor data
            df_biomass = biomass_analysis.separate_cat_ana_anycofactor2(
                model,
                macromolecule_exch_name=self.biomass_reac,
                catana_exchange_reac={k: v},
                last_cofactor=booll,
                optimize_flux=self.objective_flux
            )
        
        # Calculate ATP and cofactor-related data
        df_ATP = biomass_analysis.calculate_energy(split_reacs=self.reversiblereacs)
        df_cofactors = biomass_analysis.calculate_cofactors(split_reacs=self.reversiblereacs)
        
        return df_biomass, df_ATP, df_cofactors, biomass_analysis

    def define_required_exchanges(self, active_fluxdict, biomass_analysis, df_result_allcofactors, currencies):
        """
        Determines the required exchange reactions for the anabolic model and calculates related fluxes.
    
        This method analyzes the anabolic model to identify necessary exchange reactions based on flux variability
        and cofactor currencies. It then computes and adjusts the fluxes related to these reactions and generates 
        dataframes for biomass, ATP, and cofactors.
    
        Parameters
        ----------
        active_fluxdict : dict
            A dictionary of active fluxes where the keys are flux identifiers and the values are flux rates.
            
        biomass_analysis : macromolecule_analysis
            An analysis object used to analyze the biomass reaction and the anabolic model.
            
        df_result_allcofactors : pandas.DataFrame
            A dataframe containing cofactor normalization data for reactions.
            
        currencies : dict
            A dictionary of cofactor currencies where keys are cofactor identifiers and values are associated data.
    
        Returns
        -------
        df : pandas.DataFrame
            A dataframe with normalized flux values, including biomass, anabolic, and catabolic fluxes.
            
        atp_df : pandas.DataFrame
            A dataframe with calculated ATP-related data.
            
        cofactor_df : pandas.DataFrame
            A dataframe with calculated cofactor-related data.
            
        required_dict : dict
            A dictionary of required exchange reactions and their associated currency data.
        """
        
        # Clone the anabolic model for analysis
        anabolic_model = biomass_analysis.model_anabolic.clone()
        
        # List of reactions to analyze
        r_list = list(currencies.keys())
        
        # Perform flux variability analysis on the selected reactions
        fva = cbmpy.CBCPLEX.cplx_FluxVariabilityAnalysis(anabolic_model, selected_reactions=r_list)
        fva_dict = dict(zip(fva[1], fva[0]))
        
        required = []
        pools = []
        
        # List of reactions with flux variability data
        demands_list = list(fva_dict.keys())
        
        for k in demands_list:
            # Skip if the reaction is already in a pool
            if k in [item for sublist in pools for item in sublist]:
                continue
            
            v = fva_dict[k]
            if abs(v[3]) > 1e-6:  # Check if maximum flux is significant
                required.append(k)
                pools.append([k])
                continue
            
            # Save old bounds and set new bounds for analysis
            old_bounds = [anabolic_model.getReaction(k).getLowerBound(), anabolic_model.getReaction(k).getUpperBound()]
            anabolic_model.setReactionBounds(k, v[2], v[2])
            
            # Perform a new flux variability analysis with updated bounds
            fva_new = cbmpy.CBCPLEX.cplx_FluxVariabilityAnalysis(anabolic_model, selected_reactions=r_list)
            fva_dict_new = dict(zip(fva_new[1], fva_new[0]))
            
            # Identify reactions with no variability and flux close to zero
            specific_pool = [k]
            for kn, vn in fva_dict_new.items():
                if abs(vn[4]) < 1e-6 and abs(vn[3]) < 1e-6:
                    specific_pool.append(kn)
            
            # Restore original bounds and update pools
            anabolic_model.setReactionBounds(k, old_bounds[0], old_bounds[1])
            pools.append(specific_pool)
            required.append(specific_pool[0])
        
        # Adjust the model to include only the required reactions
        for p in pools:
            if len(p) > 1:
                for r in p[1:]:
                    anabolic_model.setReactionBounds(r, 0, 0)
        
        # Analyze the modified model
        result = cbmpy.CBGLPK.glpk_analyzeModel(anabolic_model, method='s')
        anafluxes = anabolic_model.getReactionValues()
        
        # Create a dataframe with biomass, anabolic, and catabolic fluxes
        df = pd.DataFrame.from_records([active_fluxdict, anafluxes], index=['biomass', 'anabolic']).T
        df.fillna(value=0, inplace=True)
        df['catabolic'] = np.zeros(len(df))
        
        # Compute catabolic fluxes for required reactions
        for reac in required:
            df[reac+'_norm'] = df_result_allcofactors[reac+'_norm']
            df[reac] = df[reac+'_norm'] * np.abs(anafluxes[reac]) / df.loc[reac, reac+'_norm']
            df['catabolic'] += df[reac]
        
        # Check sum consistency
        df['sumcheck'] = df['biomass'] - df['catabolic'] - df['anabolic']
        
        # Update biomass analysis with results
        biomass_analysis.df = df
        biomass_analysis.exchanged_cofactors = required
        
        # Calculate ATP and cofactor-related data
        atp_df = biomass_analysis.calculate_energy(split_reacs=self.reversiblereacs)
        cofactor_df = biomass_analysis.calculate_cofactors(split_reacs=self.reversiblereacs)
        
        # Create a dictionary of required reactions and their currency data
        required_dict = {r: currencies[r] for r in required}
        
        return df, atp_df, cofactor_df, required_dict

    def do_split_macromolecules(self, df_result, currencies):
        """
        Splits macromolecule reactions and calculates associated fluxes for the given model.
    
        This method performs an analysis of macromolecule reactions, splits them based on provided currencies, and
        calculates related fluxes, ATP, and cofactor data for each macromolecule exchange.
    
        Parameters
        ----------
        df_result : pandas.DataFrame
            A dataframe containing flux values including biomass and anabolic fluxes.
            
        currencies : dict
            A dictionary of cofactor currencies where keys are cofactor identifiers and values are associated data.
    
        Returns
        -------
        macromolecule_dfs : list of pandas.DataFrame
            A list of dataframes with separated fluxes for each macromolecule exchange.
            
        atp_macromolecule_dfs : list of pandas.DataFrame
            A list of dataframes with ATP-related data for each macromolecule exchange.
            
        cofactor_macromolecule_dfs : list of pandas.DataFrame
            A list of dataframes with cofactor-related data for each macromolecule exchange.
        """
        
        # Convert biomass fluxes to a dictionary
        fluxdict = df_result['biomass'].to_dict()
        
        macromolecule_dfs = []
        atp_macromolecule_dfs = []
        cofactor_macromolecule_dfs = []
        
        for i, exch in enumerate(self.macromolecule_exch_ids):
            # Create an analysis object for each macromolecule exchange
            analysis = macromolecule_analysis(
                self.active_model.clone(), 
                self.macromolecule_species_ids[i], 
                biomass_reac=self.biomass_reac
            )
            
            # Add reactions and species for the macromolecule exchange
            reactions = [[exch, self.macromolecule_species_ids[i] + ' exchange', {self.macromolecule_species_ids[i]: -1.0}]]
            species = []
            analysis.add_reactions_and_species(reactions, species)
            
            # Perform analysis on the macromolecule exchange reaction
            analysis.perform_analysis({exch: 1.0}, optimize_flux=self.objective_flux)
    
            # Generate separation dataframe
            df = analysis.generate_separation_df(fluxdict)
            
            # Separate reactions based on cofactors
            count = 0
            booll = False
            for k, v in currencies.items():
                count += 1
                if count >= len(currencies):
                    booll = True
                
                df = analysis.separate_cat_ana_anycofactor2(
                    active_model=analysis.active_model, 
                    biomass_anabolic=df_result['anabolic'],
                    macromolecule_exch_name=exch, 
                    catana_exchange_reac={k: v}, 
                    last_cofactor=booll, 
                    optimize_flux=self.objective_flux
                )
            
            # Calculate ATP and cofactor data for the macromolecule exchange
            atp_df = analysis.calculate_energy(split_reacs=self.reversiblereacs)
            cofactors_df = analysis.calculate_cofactors(split_reacs=self.reversiblereacs)
            
            # Append results to lists
            macromolecule_dfs.append(df)
            atp_macromolecule_dfs.append(atp_df)
            cofactor_macromolecule_dfs.append(cofactors_df)
        
        # Perform GAM analysis
        GAM_analysis = macromolecule_analysis(
            self.active_model.clone(), 
            'GAM', 
            biomass_reac=self.biomass_reac
        )
        reactions = [['R_GAM', 'Growth associated maintenance', {
            'M_atp_c': -1.0, 
            'M_h2o_c': -1.0, 
            'M_adp_c': 1.0, 
            'M_pi_c': 1.0, 
            'M_h_c': 1.0
        }]]
        species = []
        GAM_analysis.add_reactions_and_species(reactions, species)
        GAM_analysis.perform_analysis({self.GAM_id: self.GAM_value}, optimize_flux=self.objective_flux)
        
        # Generate separation dataframe for GAM
        df_GAM = GAM_analysis.generate_separation_df(fluxdict)
        count = 0
        booll = False
        for k, v in currencies.items():
            count += 1
            if count >= len(currencies):
                booll = True
            
            df_GAM = GAM_analysis.separate_cat_ana_anycofactor2(
                active_model=GAM_analysis.active_model, 
                biomass_anabolic=df_result['anabolic'],
                macromolecule_exch_name='R_GAM', 
                catana_exchange_reac={k: v}, 
                last_cofactor=booll, 
                optimize_flux=self.objective_flux
            )
        
        # Calculate ATP and cofactor data for GAM
        GAM_atp_df = GAM_analysis.calculate_energy(split_reacs=self.reversiblereacs)
        GAM_cofactors_df = GAM_analysis.calculate_cofactors(split_reacs=self.reversiblereacs)
        
        # Append GAM results to lists
        macromolecule_dfs.append(df_GAM)
        atp_macromolecule_dfs.append(GAM_atp_df)
        cofactor_macromolecule_dfs.append(GAM_cofactors_df)
        
        return macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs


    def get_activeN(self, model):
        """
        Retrieves the active stoichiometric matrix from a given genome-scale metabolic model.
        
        Parameters
        ----------
        model : cbmpy.CBModel
            The genome-scale metabolic model from which the stoichiometric matrix will be extracted.
        
        Returns
        -------
        N_active : np.array
            The active stoichiometric matrix with only active metabolites and fluxes.
        """
        
        # Get reaction values from the model
        values = model.getReactionValues()
        # Build the stoichiometric matrix from the model
        N = np.array(model.buildStoichMatrix(matrix_type='numpy', only_return=True).array)
        
        # List to store indices of inactive fluxes
        i_inactive = []
        # Identify inactive fluxes based on their values (near zero)
        for i in range(len(list(values.values()))):
            if abs(list(values.values())[i]) < 1e-10:
                i_inactive.append(i)
        
        # Remove columns corresponding to inactive fluxes from the stoichiometric matrix
        N_activeflux = np.delete(N, i_inactive, 1)
        
        # List to store indices of inactive metabolites
        i_inactive = []
        # Identify inactive metabolites (rows with all near-zero values)
        for i in range(N_activeflux.shape[0]):
            if np.all(np.abs(N_activeflux[i, :]) < 1e-10):
                i_inactive.append(i)
        
        # Remove rows corresponding to inactive metabolites from the matrix
        N_active = np.delete(N_activeflux, i_inactive, 0)
        
        # Return the active stoichiometric matrix
        return N_active
    
    def get_activeFluxdict(self, model):
        """
        Returns a dictionary of active fluxes from the model, with flux names as keys and flux values as values.
        
        Parameters
        ----------
        model : cbmpy.CBModel
            The genome-scale metabolic model from which active fluxes will be extracted.
        
        Returns
        -------
        active_fluxes : dict
            A dictionary containing the names and values of active fluxes (fluxes with values greater than 1e-10).
        """
        
        # Get the dictionary of reaction values (fluxes)
        fluxes = model.getReactionValues()
        
        # Initialize an empty dictionary to store active fluxes
        active_fluxes = {}
        # Populate the dictionary with fluxes that are significantly non-zero
        for key in fluxes:
            if abs(fluxes[key]) > 1e-10:
                active_fluxes[key] = fluxes[key]
        
        # Return the dictionary of active fluxes
        return active_fluxes
    
    def get_activemodel(self, model):
        """
        Generates a new model with only active fluxes and metabolites.
        
        Parameters
        ----------
        model : cbmpy.CBModel
            The genome-scale metabolic model from which the active model will be derived.
        
        Returns
        -------
        active_model : cbmpy.CBModel
            A new model containing only active fluxes and active metabolites.
        """
        
        # Clone the original model to create an active model
        active_model = model.clone()
        
        # Build the stoichiometric matrix from the cloned model
        N = np.array(active_model.buildStoichMatrix(matrix_type='numpy', only_return=True).array)
        # Get the reaction identifiers
        reacts = active_model.buildStoichMatrix(matrix_type='numpy', only_return=True).col
        # Get the metabolite identifiers
        metabs = active_model.buildStoichMatrix(matrix_type='numpy', only_return=True).row
        
        # Get the solution vector (flux values)
        fluxes = active_model.getSolutionVector()
        
        # List to store indices of inactive fluxes
        i_inactive = []
        # Identify inactive fluxes and delete them from the active model
        for i in range(len(fluxes)):
            if abs(fluxes[i]) < 1e-10:
                active_model.deleteReactionAndBounds(reacts[i])
                i_inactive.append(i)
        
        # Remove columns corresponding to inactive fluxes from the stoichiometric matrix
        N_activeflux = np.delete(N, i_inactive, 1)
        
        # Identify inactive metabolites and delete them from the active model
        for i in range(len(metabs)):
            if np.all(np.abs(N_activeflux[i, :]) < 1e-10):
                active_model.deleteSpecies(metabs[i])
        
        # Return the active model with only active fluxes and metabolites
        return active_model
    
    
    def get_activeconstraints(self, model, active_fluxes):
        """
        Gets the active constraints from the model based on nonzero flux values.
        
        Parameters
        ----------
        model : cbmpy.CBModel
            The genome-scale metabolic model from which active constraints will be extracted.
        active_fluxes : dict
            Dictionary with nonzero fluxes. Can be obtained using the function get_activeFluxdict.
            Fluxes with 'soft-knockout' (bounds to [0,0]) are not considered active.
        
        Returns
        -------
        active_constraints : dict
            A dictionary with reaction names and values for fluxes that hit their bounds.
        """
        
        # Lists to store names and values of constraints
        names_min = []
        Js_min = []
        
        # Loop over active fluxes
        for i in range(len(list(active_fluxes.values()))):
            # Get the flux value and its name
            J = list(active_fluxes.values())[i]
            name = list(active_fluxes.keys())[i]
            # Get the bounds of the reaction
            Upbound = model.getReactionUpperBound(name)
            Lowbound = model.getReactionLowerBound(name)
            # Check if the flux reaches either bound
            if np.abs(J - Upbound) < 1e-10 or np.abs(J - Lowbound) < 1e-10:
                names_min.append(name)
                Js_min.append(J)
          
        # Return dictionary with active constraint reaction names and flux values
        return dict(zip(names_min, Js_min))
    
    def get_rankcheck(self, model, nonzero_dict, N_active, numerical_check=True):
        """
        Prints statements about the rank, columns, the number of constraints hit, and the null space of the active stoichiometric matrix.
        Note: Rank(N) + Nconstraints = NcolumnsN only when the minimum sum of absolute fluxes has been taken.
    
        Parameters
        ----------
        model : cbmpy.CBModel
            The genome-scale metabolic model to check.
        nonzero_dict : dict
            Dictionary with nonzero fluxes. Can be obtained using the function get_activeFluxdict.
        N_active : np.array
            Array with the active stoichiometric matrix. Can be obtained using the function get_activeN.
        numerical_check : bool, optional
            Whether to perform a numerical consistency check. Default is True.
    
        Returns
        -------
        None
        """
        
        # Lists to store names and values of constraints
        names_min = []
        Js_min = []
        
        # Loop over nonzero fluxes
        for i in range(len(list(nonzero_dict.values()))):
            # Get the flux value and its name
            J = list(nonzero_dict.values())[i]
            name = list(nonzero_dict.keys())[i]
            # Get the bounds of the reaction
            Upbound = model.getReactionUpperBound(name)
            Lowbound = model.getReactionLowerBound(name)
            # Check if the flux reaches either bound
            if np.abs(J - Upbound) < 1e-10 or np.abs(J - Lowbound) < 1e-10:
                names_min.append(name)
                Js_min.append((J, name))
        
        # Calculate the rank of the active stoichiometric matrix
        rank_N_active = np.linalg.matrix_rank(N_active)
        # Calculate the null space of the active stoichiometric matrix
        null_N_active = sp.linalg.null_space(N_active)
        
        # Print the number of constraints hit
        print('Number of constraints hit: {}'.format(len(Js_min)))
        # Print the rank of the stoichiometric matrix
        print('Rank of N: {}'.format(rank_N_active))
        # Print the number of columns in the stoichiometric matrix
        print('Nr of columns in N: {}'.format(N_active.shape[1]))
        # Print the number of columns in the null space of the stoichiometric matrix
        print('Nr of columns in null space of N: {}'.format(null_N_active.shape[1]))
        
        # Perform a numerical consistency check if required
        if len(Js_min) + rank_N_active != N_active.shape[1] and numerical_check is True:
            print(Js_min)
            raise Exception('The number of active constraints plus the rank of the active stoichiometric matrix is not equal to the number of columns in the active stoichiometric matrix. Check the model for inconsistencies.')
        
        return None
    
    def generate_EFM(self, active_model, bounded_flux, bounds, unbounded_fluxes, optimize_flux, optimize_flux_bounds, sense='minimize'):
        """
        Generates an EFM from an active model (after minimization of the sum of absolute fluxes), where all inactive fluxes and metabolites are deleted.
        Expanded to accept models with more than two EFMs.
    
        Parameters
        ----------
        active_model : cbmpy.CBModel
            Model without inactive fluxes and metabolites. Ideally after minimization of the sum of absolute fluxes.
        bounded_flux : str
            Flux name of the flux to which the EFM should be related.
        bounds : list
            List with [lower_bound, upper_bound] values. Usually, lower_bound = upper_bound.
        unbounded_fluxes : list of str
            List of flux names from which the bounds should be set to [0, 1000].
        optimize_flux : str
            Flux name that should be optimized.
        optimize_flux_bounds : list
            Bounds [lower_bound, upper_bound] of the flux to be optimized.
        sense : str, optional
            Whether the flux should be minimized or maximized. Options are 'minimize' or 'maximize'. Default is 'minimize'.

        Returns
        -------
        model : cbmpy.CBModel
            The EFM model.
        """
        
        
        # Clone the active model to generate the EFM model
        model = active_model
        
        # Set bounds for the flux to be optimized
        model.setReactionBounds(optimize_flux, lower=0, upper=1000)
        model.setReactionBounds(bounded_flux, lower=bounds[0], upper=bounds[1])
        
        # Get all reaction names in the model
        unbounded_fluxes = list(model.getReactionValues().keys())
        # Remove the bounded flux from the list of unbounded fluxes
        unbounded_fluxes.remove(bounded_flux)
        # Set bounds for all unbounded fluxes
        for f in unbounded_fluxes:
            model.setReactionBounds(f, lower=0, upper=1000)
        
        # Set the objective flux and the optimization sense (minimize or maximize)
        model.setObjectiveFlux(optimize_flux, osense=sense)
        
        # Analyze the model using GLPK solver
        result = cbmpy.CBGLPK.glpk_analyzeModel(model)
        
        # Return the generated EFM model
        return model

    def get_MCEQ(self, model, norm_metab=None, norm_value=None, biomass_reac='R_BIOMASS_Ec_iML1515_core_75p37M', additional_active_metabolites=[], flux_dist=None):
        """
        Get the macrochemical equation from a model. Currently implemented with the core reaction of the E. coli model.
    
        Parameters
        ----------
        model : cbmpy.CBModel
            The genome-scale metabolic model from which the macrochemical equation will be extracted.
        norm_metab : str, optional
            Metabolite to which the macrochemical equation should be normalized. Default is None.
        norm_value : float, optional
            Value to which the macrochemical equation should be normalized. Default is None.
        extraconstraint : bool, optional
            Whether a user constraint is added in the model. In this case, cbmpy with GLPK cannot be used and CPLEX is used. Default is False.
        biomass_reac : str, optional
            The name of the biomass reaction. Default is 'R_BIOMASS_Ec_iML1515_core_75p37M'.
        additional_active_metabolites : list, optional
            List of additional active metabolites to be considered. Default is an empty list.
        flux_dist : list, optional
            List of flux distributions. Default is None.
        
        Returns
        -------
        macrochemical : np.array with sympy input
            Macrochemical equation from the model.
        """
        
        # Clone the original model to create a working model
        MCEQ_model = model.clone()
        
        # Split reversible reactions in the model
        MCEQ_model = cbmpy.CBTools.splitReversibleReactions(MCEQ_model)
        
        # Analyze the model
        result = cbmpy.CBGLPK.glpk_analyzeModel(MCEQ_model)
        
        # Get the stoichiometric matrix and reaction/metabolite identifiers
        N_min = np.array(MCEQ_model.buildStoichMatrix(matrix_type='numpy', only_return=True).array)
        reacts = MCEQ_model.buildStoichMatrix(matrix_type='numpy', only_return=True).col
        metabs = MCEQ_model.buildStoichMatrix(matrix_type='numpy', only_return=True).row
        
        # Use provided flux distribution or get the solution vector from the model
        if flux_dist is None:
            fluxes = MCEQ_model.getSolutionVector()
        else:
            fluxes = flux_dist[:len(reacts)] # if fluxes were added, delete them (such as ATP demand etc)
        
        # Identify inactive fluxes
        i_inactive_fluxes = [i for i in range(len(fluxes)) if abs(fluxes[i]) < 1e-10]
        active_fluxes_min = [reacts[i] for i in range(len(fluxes)) if abs(fluxes[i]) >= 1e-10]
    
        # Remove columns corresponding to inactive fluxes
        N_activeflux_min = np.delete(N_min, i_inactive_fluxes, 1)
    
        # Identify inactive metabolites
        i_inactive_metabs = [i for i in range(N_activeflux_min.shape[0]) if np.all(np.abs(N_activeflux_min[i, :]) < 1e-10)]
        active_metabolites_min = [metabs[i] for i in range(N_activeflux_min.shape[0]) if not np.all(np.abs(N_activeflux_min[i, :]) < 1e-10)]
        
        # Remove rows corresponding to inactive metabolites
        N_active_min = np.delete(N_activeflux_min, i_inactive_metabs, 0)
    
        # Create sympy symbols for active extracellular metabolites or additional active metabolites
        metabs_symbols = sympy.zeros(len(active_metabolites_min), 1)
        for i, metab in enumerate(active_metabolites_min):
            if metab.endswith('_e') or metab in additional_active_metabolites:
                metabs_symbols[i] = sympy.Symbol(metab)
    
        # Identify exchange reactions
        i_ex = [i for i in range(len(active_fluxes_min)) if active_fluxes_min[i].startswith('R_EX')]
        
        # Get non-zero fluxes
        fluxvec_min_notzero = [fluxes[i] for i in range(len(fluxes)) if abs(fluxes[i]) > 1e-10]
    
        # Remove columns corresponding to exchange reactions from the active stoichiometric matrix and flux vector
        N_mac = np.delete(N_active_min, i_ex, 1)
        fluxvec_mac = np.delete(fluxvec_min_notzero, i_ex)
        
        # Calculate the macrochemical equation
        macrochemical = metabs_symbols.T @ N_mac @ fluxvec_mac
        
        # Normalize the macrochemical equation if required
        if norm_value is not None and norm_metab is not None:
            metabs = []
            floats = []
            
            for arg in macrochemical[0].args:
                floats.append(np.float64(arg.args[0]))
                metabs.append(str(arg.args[1].name))
            
            ind = np.where(norm_metab == np.array(metabs))[0][0]
            macrochemical = macrochemical * norm_value / floats[ind]
        elif norm_value is not None or norm_metab is not None:
            raise Exception('Both norm_metab and norm_value must be provided')
        
        return macrochemical


class macromolecule_analysis():
    def __init__(self, active_model_org, molecule_name, biomass_reac):
        """
        Initialize the MacromoleculeAnalysis class.
        
        Parameters:
        active_model_org : cbmpy.CBModel
            The original active model to be cloned.
        molecule_name : str
            The name of the molecule.
        biomass_reac : str
            The biomass reaction ID.
        """
        self.name = molecule_name
        self.active_model = active_model_org.clone()  # Clone the original active model
        self.objective_value = None
        self.active_fluxdict = None
        self.MCEQ = None
        self.MCEQ_selected = None
        self.N_active = None
        self.active_constraints = None
        self.df = None
        self.energy_obj = None
        self.energy_df = None
        self.cofactor_df = None
        self.model_anabolic = None
        self.exchanged_cofactors = []
        self.exch_reacs = []
        self.biomass_reac = biomass_reac  # Biomass reaction ID
        self.active_model.setReactionBounds(self.biomass_reac, 0, 1000)  # Set biomass reaction bounds

    def add_reactions_and_species(self, reactions, species):
        """
        Add reactions and species to the model.
        
        Parameters:
        reactions : list of lists
            Each list contains [id, name, metabolite:coeff dict].
        species : list of lists
            Each list contains [id, name, compartment].
        """
        active_model = self.active_model.clone()  # Clone the active model
        
        for i in species:
            active_model.createSpecies(i[0], name=i[1], compartment=i[2])  # Add species to the model
        
        for i in reactions:
            active_model.createReaction(i[0], name=i[1], reversible=False)  # Add reactions to the model
            for j in list(i[2].keys()):
                active_model.createReactionReagent(i[0], metabolite=j, coefficient=i[2][j])  # Add metabolites to reactions
            active_model.setReactionBounds(i[0], 0, 1000)  # Set reaction bounds
        
        self.active_model = active_model.clone()  # Update the active model

    def perform_analysis(self, flux_bounds_lb_spec, optimize_flux='R_EX_glc__D_e_rev'):
        """
        Perform analysis on the active model.
        
        Parameters:
        flux_bounds_lb_spec : dict
            Lower bounds for flux constraints.
        optimize_flux : str
            The flux to optimize.
        """
        ama = active_model_analysis(self.active_model)  # Initialize active model analysis
        ama.perform_analysis(flux_bounds_lb=flux_bounds_lb_spec, optimize_flux=optimize_flux)  # Perform the analysis
        self.MCEQ = ama.MCEQ  # Store MCEQ results
        self.objective_value = ama.objective_value  # Store objective value
        self.active_fluxdict = ama.active_fluxdict  # Store active flux dictionary
        self.N_active = ama.N_active  # Store active N values
        self.active_constraints = ama.active_constraints  # Store active constraints

    def generate_separation_df(self, fluxdict_biomass):
        """
        Generate a separation DataFrame.
        
        Parameters:
        fluxdict_biomass : dict
            Flux dictionary for biomass.
        
        Returns:
        DataFrame
            A DataFrame with separation values.
        """
        norm_reac = self.biomass_reac  # Normalize reaction by biomass reaction
        df = pd.DataFrame.from_records([fluxdict_biomass, self.active_fluxdict], index=['biomass', '{}+biomass'.format(self.name)]).T
        df.fillna(value=0, inplace=True)  # Fill NaN values with 0
        df[self.name] = df['{}+biomass'.format(self.name)] - df['{}+biomass'.format(self.name)].loc[norm_reac]/df['biomass'].loc[norm_reac]*df['biomass']
        self.df = df  # Store the DataFrame
        
        return df
    
    def separate_cat_ana_anycofactor2(self, active_model, biomass_anabolic=None, optimize_flux='R_EX_glc__D_e_rev', macromolecule_exch_name='R_EX_PROTEIN', 
                                      catana_exchange_reac={'R_ATP_demand':  [('M_atp_c', -1.0), ('M_h2o_c', -1.0), ('M_adp_c', 1.0), ('M_h_c', 1.0), ('M_pi_c', 1.0)]}, 
                                      last_cofactor=False):
        """
        Separate catabolic and anabolic reactions including any cofactor.
        
        Parameters:
        active_model : cbmpy.CBModel
            The active model to use.
        biomass_anabolic : dict, optional
            Anabolic biomass reactions, default is None.
        optimize_flux : str, optional
            Flux to optimize, default is 'R_EX_glc__D_e_rev'.
        macromolecule_exch_name : str, optional
            Name of the macromolecule exchange reaction, default is 'R_EX_PROTEIN'.
        catana_exchange_reac : dict, optional
            Exchange reactions, default is {'R_ATP_demand':  [('M_atp_c', -1.0), ('M_h2o_c', -1.0), ('M_adp_c', 1.0), ('M_h_c', 1.0), ('M_pi_c', 1.0)]}.
        last_cofactor : bool, optional
            Whether this is the last cofactor, default is False.
        
        Returns:
        DataFrame
            A DataFrame with separated reactions.
        """
        df = self.df  # Reference the DataFrame
        
        if self.name != 'biomass':
            if biomass_anabolic is None:
                raise Exception('If you are working with a macromolecule you need to specify the anabolic biomass reaction!')
            else:
                df['anabolic_biomass'] = df.index.map(biomass_anabolic)  # Map anabolic biomass reactions
        
        active_model_catabolic = active_model.clone()  # Clone the active model
        active_fluxdict = active_model.getReactionValues()  # Get reaction values
        
        exch_r = list(catana_exchange_reac.keys())[0]  # Get the first exchange reaction
        if exch_r not in df.index:
            df.loc[exch_r] = np.zeros(len(df.columns))  # Initialize the exchange reaction in the DataFrame
        r_descr = list(catana_exchange_reac.values())[0]  # Get the reaction description
        
        colname = exch_r

        active_model_catabolic.createReaction(exch_r, reversible=False, name='pseudoreaction for '+exch_r)  # Create exchange reaction
        
        for tup in r_descr:
            active_model_catabolic.createReactionReagent(exch_r, tup[0], tup[1])  # Add reagents to the reaction
        active_model_catabolic.setReactionBounds(exch_r, 1, 1000)  # Set reaction bounds
        
        active_model_catabolic.setReactionBounds(self.biomass_reac, 0, 1000)  # Set biomass reaction bounds
        active_model_catabolic.setReactionBounds(macromolecule_exch_name, 0, 1000)  # Set macromolecule exchange reaction bounds
        
        active_model_catabolic.setObjectiveFlux(optimize_flux, osense='minimize')  # Set the objective flux
        
        active_result_biomass_catabolic = cbmpy.CBGLPK.glpk_analyzeModel(active_model_catabolic, method='s')  # Analyze the model
        active_fluxes_biomass_catabolic = active_model_catabolic.getReactionValues()  # Get the reaction values
        
        self.exchanged_cofactors.append(colname)  # Append the cofactor
        self.exch_reacs.append(exch_r)  # Append the exchange reaction
        
        if self.name == 'biomass':
            df[colname+'+biomass_norm'] = df.index.map(active_fluxes_biomass_catabolic)  # Normalize the biomass column
            df.fillna(value=0, inplace=True)  # Fill NaN values with 0
            df[colname+'_norm'] = df[colname+'+biomass_norm'] - df[colname+'+biomass_norm'].loc[self.biomass_reac]/df['biomass'].loc[self.biomass_reac]*df['biomass']
        else:
            df[colname+'+macromol+biomass_norm'] = df.index.map(active_fluxes_biomass_catabolic)  # Normalize the macromolecule column
            df.fillna(value=0, inplace=True)  # Fill NaN values with 0
            df[colname+'+macromol_norm'] = df[colname+'+macromol+biomass_norm'] - df[colname+'+macromol+biomass_norm'].loc[self.biomass_reac]/df['biomass'].loc[self.biomass_reac]*df['biomass']
            df[colname+'_norm'] = df[colname+'+macromol_norm'] - df[colname+'+macromol_norm'].loc[macromolecule_exch_name]/df[self.name].loc[macromolecule_exch_name]*df[self.name]
        
        df[colname] = df.index.map(active_fluxes_biomass_catabolic)  # Map the active fluxes
        
        #TODO test this and the version below (from the old version, should be commented. I don't know why it's so different)
        # active_model_anabolic = active_model.clone()  # Clone the active model again
        
        # for r, v in active_fluxdict.items():
        #     if r != optimize_flux:
        #         if r in df[df[self.name] < 0].index:  # Check if the reaction is anabolic
        #             active_model_anabolic.setReactionBounds(r, 0, 0)  # Set reaction bounds to zero
        #         elif 'EX_' in r and r != 'R_EX' and r != 'R_EX_' + self.name:
        #             if r not in self.exch_reacs:
        #                 active_model_anabolic.setReactionBounds(r, -1000, 0)  # Allow uptake reactions
        #             active_model_anabolic.setReactionBounds(r, 0, 0)  # Set bounds for exchange reactions
        #         else:
        #             active_model_anabolic.setReactionBounds(r, -1000, 0)  # Set bounds for other reactions
        #             active_model_anabolic.setReactionBounds(r, 0, 0)  # Reset the bounds
        #     else:
        #         active_model_anabolic.setReactionBounds(r, 0, 1000)  # Allow production reactions
        
        # active_model_anabolic.setReactionBounds(optimize_flux, -1000, 0)  # Set the objective flux bounds
        
        # if last_cofactor:
        #     self.model_anabolic = active_model_anabolic.clone()  # Store the anabolic model
        # else:
        #     self.model_anabolic = active_model_anabolic.clone()  # Clone the anabolic model
        #     self.active_fluxdict = active_model_anabolic.getReactionValues()  # Get the reaction values
        
        # return df
        
        
        
        if len(self.exchanged_cofactors) == 1:
            active_model_anabolic = active_model_catabolic.clone()
            
        else:
            
            active_model_anabolic = self.model_anabolic.clone()
            
            exch_r = list(catana_exchange_reac.keys())[0]
            r_descr = list(catana_exchange_reac.values())[0]
            colname=exch_r
            if exch_r not in active_model_anabolic.getReactionIds():
                active_model_anabolic.createReaction(exch_r, reversible=False, name='pseudoreaction for '+exch_r)
                for tup in r_descr:    
                    active_model_anabolic.createReactionReagent(exch_r, tup[0], tup[1])
            active_model_anabolic.setReactionBounds(exch_r, -1000, 0)
            
        
        reac_list = active_model_anabolic.buildStoichMatrix(matrix_type='numpy', only_return=True).col
        
        for r in reac_list:
            # print(r)
            if df.loc[r, colname+'_norm']>1e-9:
                # print(r)
                # print(colname+'_norm')
                stoich = active_model_anabolic.getReaction(r).getStoichiometry()
                chem = active_model_anabolic.getSpecies(stoich[0][1]).getChemFormula()
                if r[:4] == 'R_EX':
                    # active_model_anabolic.setReactionBounds(r, -1000, 1000)
                    if chem is None or 'C' in chem:
                        active_model_anabolic.setReactionBounds(r, 0, 1000)
                    elif r == active_model_anabolic.getActiveObjectiveReactionIds()[0]:
                        active_model_anabolic.setReactionBounds(r, 0, 1000)
                    elif r in ['R_EX_lac__D_e', 'R_EX_ac_e', 'R_EX_succ_e', 'R_EX_co2_e_fwd', 'R_EX_pyr_e']:#TODO remove later (Look at this? test with and without)
                        active_model_anabolic.setReactionBounds(r, 0, 1000)
                        # raise Exception('test')
                    else:
                        active_model_anabolic.setReactionBounds(r, -1000, 1000)
                elif r == exch_r:
                    active_model_anabolic.setReactionBounds(r, -1000, 0)
                elif (1.0, 'M_co2_c') in stoich or (1.0, 'M_co2_m') in stoich:# or 'C' in chem:
                    print(chem)
                    active_model_anabolic.setReactionBounds(r, 0, 1000)
                
                else:
                    active_model_anabolic.setReactionBounds(r, -1000, 1000)
        
        

        
        active_model_anabolic.setObjectiveFlux(optimize_flux, osense='minimize')
        active_model_anabolic.setReactionBounds(macromolecule_exch_name, df[self.name].loc[macromolecule_exch_name], df[self.name].loc[macromolecule_exch_name])
        
        active_result_biomass_anabolic = cbmpy.CBGLPK.glpk_analyzeModel(active_model_anabolic, method='s')
        
        active_fluxes_biomass_anabolic = active_model_anabolic.getReactionValues()

        
        self.model_anabolic = active_model_anabolic.clone()

        # print('ATP demand from FBA' + str(active_fluxes_biomass_anabolic['R_ATP_demand']))
        if self.name == 'biomass':
            df['anabolic'] = df.index.map(active_fluxes_biomass_anabolic)
            
            df[colname] = -1*df['anabolic'].loc[exch_r]/self.df[colname+'_norm'].loc[exch_r]*self.df[colname+'_norm']
            
            df.fillna(value=0, inplace=True)
            # df['sumcheck'] = df[self.name]-df['catabolic']-df['anabolic']
        else:
            df['anabolic_mac+biomass'] = df.index.map(active_fluxes_biomass_anabolic)
            
            df['anabolic'] = df['anabolic_mac+biomass']-df['anabolic_mac+biomass'].loc[self.biomass_reac]/df['anabolic_biomass'].loc[self.biomass_reac]*df['anabolic_biomass']
            
            df[colname] = -1*df['anabolic'].loc[exch_r]/self.df[colname+'_norm'].loc[exch_r]*self.df[colname+'_norm']
            
            
            df.fillna(value=0, inplace=True)
            
            # df['sumcheck'] = df[self.name]-df['catabolic']-df['anabolic']
        
        # print('--------------------------------------------')
        # print(self.exchanged_cofactors)
        if last_cofactor is True:
            for i, cofactor in enumerate(self.exchanged_cofactors):
                reac = self.exch_reacs[i]
                
                df[cofactor] = -1*df['anabolic'].loc[reac]/self.df[cofactor+'_norm'].loc[reac]*self.df[cofactor+'_norm']
            
            
            df['catabolic'] = np.zeros(len(df.index))
            
            for i, cofactor in enumerate(self.exchanged_cofactors):
                df['catabolic'] += df[cofactor]
            df.fillna(value=0, inplace=True)
            df['sumcheck'] = df[self.name]-df['catabolic']-df['anabolic']
        
        
        self.df=df

        return df
        
    
    def calculate_energy(self, split_reacs=True):
        """
        Calculate energy consumption and production in different metabolic pathways.
    
        Parameters
        ----------
        split_reacs : bool, optional
            Whether to split reactions into their component parts (default is True).
    
        Returns
        -------
        pd.DataFrame
            DataFrame with energy production, consumption, net energy, and special energy consumption for each pathway.
        """
        energy_obj = energy_generation(self.active_model)
        self.energy_obj = energy_obj
    
        # Initialize indices for energy DataFrame
        indices = ['tot', 'cat', 'ana'] + self.exchanged_cofactors
    
        # Create DataFrame to store energy values
        energy_df = pd.DataFrame(index=indices, columns=['prod', 'cons', 'net'])
    
        # Calculate energy values for total, catabolic, and anabolic pathways
        tot_consumption, tot_production = energy_obj.calculate_energy(self.df[self.name], split_names=split_reacs)
        cat_consumption, cat_production = energy_obj.calculate_energy(self.df['catabolic'], split_names=split_reacs)
        ana_consumption, ana_production = energy_obj.calculate_energy(self.df['anabolic'], split_names=split_reacs)
        #TODO remove this special consumption and transhydrogenase things, but must be done in the calculate_energy function etc
        # Store energy values in DataFrame
        energy_df.loc['tot', ['prod', 'cons', 'net']] = [tot_production, tot_consumption, tot_production - tot_consumption]
        energy_df.loc['cat', ['prod', 'cons', 'net']] = [cat_production, cat_consumption, cat_production - cat_consumption]
        energy_df.loc['ana', ['prod', 'cons', 'net']] = [ana_production, ana_consumption, ana_production - ana_consumption]
    
        # Calculate energy values for each exchanged cofactor
        for cofactor in self.exchanged_cofactors:
            cons, prod = energy_obj.calculate_energy(self.df[cofactor], split_names=split_reacs)
            energy_df.loc[cofactor, ['prod', 'cons', 'net']] = [prod, cons, prod - cons]
    
        self.energy_df = energy_df
        return energy_df
    
    def calculate_cofactors(self, split_reacs=True):
        """
        Calculate the production, consumption, and net balance of NADH and NADPH cofactors, including transhydrogenase activities, for various metabolic pathways and exchanged cofactors.
    
        Parameters
        ----------
        split_reacs : bool, optional
            Whether to split reactions into their component parts (default is True).
    
        Returns
        -------
        pd.DataFrame
            DataFrame with the following columns for NADH and NADPH:
            - 'prod': Production
            - 'cons': Consumption
            - 'net': Net production (prod - cons)
            - 'transhydrogenase prod': Transhydrogenase production
            - 'transhydrogenase cons': Transhydrogenase consumption
        """
        energy_obj = self.energy_obj
    
        # Define the indices for the DataFrame
        indices = ['tot', 'cat', 'ana'] + self.exchanged_cofactors
        cofactor_df = pd.DataFrame(index=indices, columns=[
            'NADH prod', 'NADH cons', 'NADH net',
            'NADPH prod', 'NADPH cons', 'NADPH net'
        ])
    
        # Calculate NADH production and consumption for total, catabolic, and anabolic pathways
        tot_nadh_prod, tot_nadh_cons = energy_obj.calculate_NADH_prod(self.df[self.name], split_names=split_reacs)
        cat_nadh_prod, cat_nadh_cons = energy_obj.calculate_NADH_prod(self.df['catabolic'], split_names=split_reacs)
        ana_nadh_prod, ana_nadh_cons = energy_obj.calculate_NADH_prod(self.df['anabolic'], split_names=split_reacs)
    
        # Store NADH values in DataFrame
        cofactor_df.loc['tot', 'NADH prod'] = tot_nadh_prod
        cofactor_df.loc['cat', 'NADH prod'] = cat_nadh_prod
        cofactor_df.loc['ana', 'NADH prod'] = ana_nadh_prod
    
        cofactor_df.loc['tot', 'NADH cons'] = tot_nadh_cons
        cofactor_df.loc['cat', 'NADH cons'] = cat_nadh_cons
        cofactor_df.loc['ana', 'NADH cons'] = ana_nadh_cons
    
        cofactor_df.loc['tot', 'NADH net'] = tot_nadh_prod - tot_nadh_cons
        cofactor_df.loc['cat', 'NADH net'] = cat_nadh_prod - cat_nadh_cons
        cofactor_df.loc['ana', 'NADH net'] = ana_nadh_prod - ana_nadh_cons
    
        # Calculate NADH values for exchanged cofactors
        for cofactor in self.exchanged_cofactors:
            nadh_prod, nadh_cons = energy_obj.calculate_NADH_prod(self.df[cofactor], split_names=split_reacs)
            cofactor_df.loc[cofactor, 'NADH prod'] = nadh_prod
            cofactor_df.loc[cofactor, 'NADH cons'] = nadh_cons
            cofactor_df.loc[cofactor, 'NADH net'] = nadh_prod - nadh_cons
    
        # Calculate NADPH production and consumption for total, catabolic, and anabolic pathways
        tot_nadph_prod, tot_nadph_cons = energy_obj.calculate_NADPH_prod(self.df[self.name], split_names=split_reacs)
        cat_nadph_prod, cat_nadph_cons = energy_obj.calculate_NADPH_prod(self.df['catabolic'], split_names=split_reacs)
        ana_nadph_prod, ana_nadph_cons = energy_obj.calculate_NADPH_prod(self.df['anabolic'], split_names=split_reacs)
    
        # Store NADPH values in DataFrame
        cofactor_df.loc['tot', 'NADPH prod'] = tot_nadph_prod
        cofactor_df.loc['cat', 'NADPH prod'] = cat_nadph_prod
        cofactor_df.loc['ana', 'NADPH prod'] = ana_nadph_prod
    
        cofactor_df.loc['tot', 'NADPH cons'] = tot_nadph_cons
        cofactor_df.loc['cat', 'NADPH cons'] = cat_nadph_cons
        cofactor_df.loc['ana', 'NADPH cons'] = ana_nadph_cons
    
        cofactor_df.loc['tot', 'NADPH net'] = tot_nadph_prod - tot_nadph_cons
        cofactor_df.loc['cat', 'NADPH net'] = cat_nadph_prod - cat_nadph_cons
        cofactor_df.loc['ana', 'NADPH net'] = ana_nadph_prod - ana_nadph_cons
    
        # Calculate NADPH values for exchanged cofactors
        for cofactor in self.exchanged_cofactors:
            nadph_prod, nadph_cons = energy_obj.calculate_NADPH_prod(self.df[cofactor], split_names=split_reacs)
            cofactor_df.loc[cofactor, 'NADPH prod'] = nadph_prod
            cofactor_df.loc[cofactor, 'NADPH cons'] = nadph_cons
            cofactor_df.loc[cofactor, 'NADPH net'] = nadph_prod - nadph_cons
           
        self.cofactor_df = cofactor_df
        return cofactor_df

    # def get_MCEQ_selected(self, fluxvec, norm_metab=None, norm_value=None):
    #     """
    #     Get the macrochemical equation from a model, specifically for the core reaction of the E. coli model.
    
    #     Parameters
    #     ----------
    #     fluxvec : pd.Series
    #         Series containing the flux values for the model.
    #     norm_metab : str, optional
    #         Metabolite to which the macrochemical equation should be normalized.
    #     norm_value : float, optional
    #         Value to which the macrochemical equation should be normalized.
    
    #     Returns
    #     -------
    #     sympy.Matrix
    #         Macrochemical equation from the model, in sympy format.
        
    #     Raises
    #     ------
    #     Exception
    #         If both norm_metab and norm_value are not provided or are mismatched.
    #     """
    #     MCEQ_model = self.active_model
        
    #     # Split reversible reactions into two unidirectional reactions
    #     MCEQ_model = cbmpy.CBTools.splitReversibleReactions(MCEQ_model)
    
    #     # Analyze the model and obtain the stoichiometric matrix
    #     result = cbmpy.CBGLPK.glpk_analyzeModel(MCEQ_model)
    #     N_min = np.array(MCEQ_model.buildStoichMatrix(matrix_type='numpy', only_return=True).array)
    #     reacts = MCEQ_model.buildStoichMatrix(matrix_type='numpy', only_return=True).col
    #     metabs = MCEQ_model.buildStoichMatrix(matrix_type='numpy', only_return=True).row
    
    #     # Filter out inactive fluxes
    #     i_inactive_fluxes = []
    #     active_fluxes_min = []
    #     to_delete = []
    
    #     for i, r in enumerate(fluxvec.index):
    #         if r not in reacts:
    #             to_delete.append(i)
        
    #     fluxvec = np.array(fluxvec)
    #     fluxvec = np.delete(fluxvec, to_delete)
    
    #     for i in range(len(fluxvec)):
    #         if abs(fluxvec[i]) < 1e-10:
    #             i_inactive_fluxes.append(i)
    #         else:
    #             active_fluxes_min.append(reacts[i])
    
    #     N_activeflux_min = np.delete(N_min, i_inactive_fluxes, 1)
    
    #     # Filter out inactive metabolites
    #     i_inactive_metabs = []
    #     active_metabolites_min = []
    
    #     for i in range(N_activeflux_min.shape[0]):
    #         if np.all(np.abs(N_activeflux_min[i, :]) < 1e-10):
    #             i_inactive_metabs.append(i)
    #         else:
    #             active_metabolites_min.append(metabs[i])
    
    #     N_active_min = np.delete(N_activeflux_min, i_inactive_metabs, 0)
    
    #     # Create symbolic representation for active metabolites
    #     metabs_symbols = sympy.zeros(len(active_metabolites_min), 1)
    #     for i in range(len(active_metabolites_min)):
    #         if active_metabolites_min[i][-2:] == '_e':
    #             metabs_symbols[i] = sympy.Symbol(active_metabolites_min[i])
    
    #     # Prepare flux vector and stoichiometric matrix for macrochemical equation
    #     fluxvec_min_notzero = []
    #     i_zero = []
    
    #     for i in range(len(fluxvec)):
    #         if abs(fluxvec[i]) > 1e-10:
    #             fluxvec_min_notzero.append(fluxvec[i])
    #             i_zero.append(i)
    
    #     interesting_reacs = ['R_EX_glc__D_e_rev', 'R_EX_BIOMASS', 'R_EX_co2_e_fwd', 'R_EX_h2o_e_fwd', 
    #                          'R_EX_h_e_fwd', 'R_EX_nh4_e_rev', 'R_EX_o2_e_rev', 'R_EX_pi_e_rev', 'R_EX_PROTEIN',
    #                          'R_EX_LIPID', 'R_EX_LPS', 'R_EX_MUREIN', 'R_EX_DNA', 'R_EX_RNA', 'R_EX_ION', 'R_GAM',
    #                          'R_EX_COFACTOR', 'R_EX_ac_e_rev', 'R_EX_etoh_e_rev', 'R_EX_for_e_rev', 'R_EX_succ_e_rev',
    #                          'R_EX_ac_e_fwd', 'R_EX_etoh_e_fwd', 'R_EX_for_e_fwd', 'R_EX_succ_e_fwd', 'R_EX_so4_e_rev']
        
    #     i_ex = []
    #     for i in range(len(active_fluxes_min)):
    #         if active_fluxes_min[i] in interesting_reacs:
    #             i_ex.append(i)
    
    #     # Compute the macrochemical equation
    #     N_mac = np.delete(N_active_min, i_ex, 1)
    #     fluxvec_mac = np.delete(fluxvec_min_notzero, i_ex)
    #     macrochemical = metabs_symbols.T @ N_mac @ fluxvec_mac
    
    #     # Normalize the macrochemical equation if required
    #     if norm_value is not None and norm_metab is not None:
    #         metabs = []
    #         floats = []
    #         for i in range(len(macrochemical[0].args)):
    #             floats.append(np.float64(macrochemical[0].args[i].args[0]))
    #             metabs.append(str(macrochemical[0].args[i].args[1].name))
    
    #         ind = np.where(norm_metab == np.array(metabs))[0][0]
    #         macrochemical = macrochemical * norm_value / floats[ind]
    #     elif norm_value is not None or norm_metab is not None:
    #         raise Exception('Fill in both metabolite name and value')
    
    #     # Format the macrochemical equation for better readability
    #     for i in sympy.preorder_traversal(macrochemical):
    #         if isinstance(i, sympy.Float):
    #             macrochemical[0] = macrochemical[0].subs(i, '{:0.3e}'.format(i))
    
    #     self.MCEQ_selected = macrochemical
    #     return macrochemical
#TODO remove this if it is not required

class active_model_analysis:
    """
    A class for performing and analyzing metabolic model optimizations and constraints.

    Attributes
    ----------
    active_model : cbmpy.CBModel
        The metabolic model to be analyzed.
    objective_value : float, optional
        The value of the objective function after optimization.
    active_fluxdict : dict, optional
        Dictionary of active fluxes and their values.
    MCEQ : sympy.Matrix, optional
        Macrochemical equation of the model.
    N_active : np.array, optional
        Active stoichiometric matrix of the model.
    active_constraints : dict, optional
        Dictionary of active constraints in the model.
    """

    def __init__(self, active_model):
        """
        Initializes the active_model_analysis with a metabolic model.

        Parameters
        ----------
        active_model : cbmpy.CBModel
            The metabolic model to be analyzed.
        """
        self.active_model = active_model
        self.objective_value = None
        self.active_fluxdict = None
        self.MCEQ = None
        self.N_active = None
        self.active_constraints = None

    def perform_analysis(self, optimize_flux='R_EX_glc__D_e_rev', sense='minimize', flux_bounds_lb={'R_EX_BIOMASS':1.0}, flux_bounds_ub={}):
        """
        Performs optimization analysis on the metabolic model and calculates active fluxes, constraints, and the macrochemical equation.

        Parameters
        ----------
        optimize_flux : str, optional
            The reaction to be optimized (default is 'R_EX_glc__D_e_rev').
        sense : str, optional
            Optimization sense ('minimize' or 'maximize', default is 'minimize').
        flux_bounds_lb : dict, optional
            Dictionary of lower bounds for fluxes (default is {'R_EX_BIOMASS': 1.0}).
        flux_bounds_ub : dict, optional
            Dictionary of upper bounds for fluxes (default is empty).

        Returns
        -------
        None
        """
        # Set lower bounds for fluxes
        for flux, lb in flux_bounds_lb.items():
            self.active_model.setReactionBounds(flux, lb, 1000)
        
        # Set upper bounds for fluxes
        for flux, ub in flux_bounds_ub.items():
            self.active_model.setReactionBounds(flux, 0, ub)
        
        # Set the objective function for optimization
        self.active_model.setObjectiveFlux(optimize_flux, osense=sense)
        
        # Analyze the model using GLPK solver
        result = cbmpy.CBGLPK.glpk_analyzeModel(self.active_model, method='s', quiet=True)
        self.objective_value = result
        
        # Obtain active stoichiometric matrix and flux dictionary
        self.N_active = separate_cat_ana.get_activeN(self, self.active_model)
        self.active_fluxdict = separate_cat_ana.get_activeFluxdict(self, self.active_model)

        # Get rank check details and constraints
        self.get_rankcheck()
        self.active_constraints = separate_cat_ana.get_activeconstraints(self, self.active_model, self.active_fluxdict)
        
        # Obtain macrochemical equation
        self.MCEQ = separate_cat_ana.get_MCEQ(self, self.active_model)

        print(f'Analysis performed; objective value = {self.objective_value}')

    def get_rankcheck(self):
        """
        Prints statements about the rank, columns, number of constraints hit, and the null space of the active stoichiometric matrix.

        Checks the consistency of the rank of the stoichiometric matrix and the constraints.

        Returns
        -------
        None
        """
        # Lists to store names of constraints and their values
        names_min = []
        Js_min = []
        
        # Loop through the active flux dictionary
        for i in range(len(self.active_fluxdict)):
            flux_value = list(self.active_fluxdict.values())[i]
            flux_name = list(self.active_fluxdict.keys())[i]
            
            # Retrieve the bounds of the reaction
            Upbound = self.active_model.getReactionUpperBound(flux_name)
            Lowbound = self.active_model.getReactionLowerBound(flux_name)
            
            # Check if the flux is at its bounds
            if np.abs(flux_value - Upbound) < 1e-10 or np.abs(flux_value - Lowbound) < 1e-10:
                names_min.append(flux_name)
                Js_min.append(flux_value)
        
        # Calculate the rank and the null space of the active stoichiometric matrix
        rank_N_active = np.linalg.matrix_rank(self.N_active)
        null_N_active = sp.linalg.null_space(self.N_active)
        
        # Print rank, columns, and active constraints statistics
        print(f'Number of constraints hit: {len(Js_min)}')
        print(f'Rank of N: {rank_N_active}')
        print(f'Number of columns in N: {self.N_active.shape[1]}')
        print(f'Number of columns in null space of N: {null_N_active.shape[1]}')
        
        return None

class energy_generation():
    def __init__(self, model):
        """
        Initializes the energy_generation class.

        Parameters
        ----------
        model : cbmpy.CBModel
            The metabolic model to be analyzed.
        """
        self.model = model
        
        # Generate a list of reactions associated with ATP, GTP, CTP, and UTP
        self.energy_reacs = self.generate_ATPreaclist()
        
        # Generate lists of NADH and NADPH associated reactions
        self.NADH_list = self.generate_NADH_list()

        self.NADPH_list = self.generate_NADPH_list()

    def generate_ATPreaclist(self):
        """
        Generates a list of reactions in the model that are associated with ATP, GTP, CTP, and UTP.

        Returns
        -------
        energy_reacs : list
            List of unique reaction IDs that involve ATP, GTP, CTP, or UTP.
        """
        # Fetch reactions associated with each energy molecule
        atp_reacs = self.model.getReactionIdsAssociatedWithSpecies('M_atp_c')
        gtp_reacs = self.model.getReactionIdsAssociatedWithSpecies('M_gtp_c')
        ctp_reacs = self.model.getReactionIdsAssociatedWithSpecies('M_ctp_c')
        utp_reacs = self.model.getReactionIdsAssociatedWithSpecies('M_utp_c')
        #TODO add mito reactions and check if this affects results (and how?)
        # Combine all reaction IDs and remove duplicates
        energy_reacs = set()
        energy_reacs.update([r[0] for r in atp_reacs])
        energy_reacs.update([r[0] for r in gtp_reacs])
        energy_reacs.update([r[0] for r in ctp_reacs])
        energy_reacs.update([r[0] for r in utp_reacs])
        
        return list(energy_reacs)

    def calculate_energy(self, active_fluxdict, growth_condition='aerobic_glucose', split_names=True):
        """
        Calculates the energy consumption and production from active reactions based on the flux dictionary.

        Parameters
        ----------
        active_fluxdict : dict
            Dictionary with active flux names and their values.
        growth_condition : string, optional
            Specifies the growth condition. Options are 'aerobic_glucose', 'anaerobic_glucose', or 'acetate'.
            Default is 'aerobic_glucose'.
        split_names : bool, optional
            Whether to use split names for reactions. Default is True.

        Returns
        -------
        list of floats
            [energy_consumption, energy_production, special_consumption].
        """

        total_production = 0
        total_consumption = 0
        
        for r in self.energy_reacs:
            if r in active_fluxdict:
                # Get species IDs for the reaction
                s = self.model.getReaction(r).getSpeciesIds()

                # Define species based on reaction direction and split_names flag
                if r[-4:] == '_fwd' and split_names:
                    st_atp = r + '_M_atp_c_fwd'
                    st_gtp = r + '_M_gtp_c_fwd'
                    st_ctp = r + '_M_ctp_c_fwd'
                    st_utp = r + '_M_utp_c_fwd'
                    
                    st_adp = r + '_M_adp_c_fwd'
                    st_gdp = r + '_M_gdp_c_fwd'
                    st_cdp = r + '_M_cdp_c_fwd'
                    st_udp = r + '_M_udp_c_fwd'
                    
                    st_amp = r + '_M_amp_c_fwd'
                    st_gmp = r + '_M_gmp_c_fwd'
                    st_cmp = r + '_M_cmp_c_fwd'
                    st_ump = r + '_M_ump_c_fwd'
                elif r[-4:] == '_rev' and split_names:
                    st_atp = r + '_M_atp_c_rev'
                    st_gtp = r + '_M_gtp_c_rev'
                    st_ctp = r + '_M_ctp_c_rev'
                    st_utp = r + '_M_utp_c_rev'
                    
                    st_adp = r + '_M_adp_c_rev'
                    st_gdp = r + '_M_gdp_c_rev'
                    st_cdp = r + '_M_cdp_c_rev'
                    st_udp = r + '_M_udp_c_rev'
                    
                    st_amp = r + '_M_amp_c_rev'
                    st_gmp = r + '_M_gmp_c_rev'
                    st_cmp = r + '_M_cmp_c_rev'
                    st_ump = r + '_M_ump_c_rev'
                else:
                    st_atp = r + '_M_atp_c'
                    st_gtp = r + '_M_gtp_c'
                    st_ctp = r + '_M_ctp_c'
                    st_utp = r + '_M_utp_c'
                    
                    st_adp = r + '_M_adp_c'
                    st_gdp = r + '_M_gdp_c'
                    st_cdp = r + '_M_cdp_c'
                    st_udp = r + '_M_udp_c'
                    
                    st_amp = r + '_M_amp_c'
                    st_gmp = r + '_M_gmp_c'
                    st_cmp = r + '_M_cmp_c'
                    st_ump = r + '_M_ump_c'
        
                c = 0
                if 'M_atp_c' in s:
                    # Calculate ATP consumption/production
                    c_atp = self.model.getReaction(r).getReagent(st_atp).getCoefficient()
                    c_axp = self.model.getReaction(r).getReagent(st_adp).getCoefficient() if 'M_adp_c' in s else \
                            self.model.getReaction(r).getReagent(st_amp).getCoefficient() if 'M_amp_c' in s else 1
                    c += c_atp

                if 'M_gtp_c' in s:
                    # Calculate GTP consumption/production
                    c_gtp = self.model.getReaction(r).getReagent(st_gtp).getCoefficient()
                    c_gxp = self.model.getReaction(r).getReagent(st_gdp).getCoefficient() if 'M_gdp_c' in s else \
                            self.model.getReaction(r).getReagent(st_gmp).getCoefficient() if 'M_gmp_c' in s else 1
                    c += c_gtp

                if 'M_utp_c' in s:
                    # Calculate UTP consumption/production
                    c_utp = self.model.getReaction(r).getReagent(st_utp).getCoefficient()
                    c_uxp = self.model.getReaction(r).getReagent(st_udp).getCoefficient() if 'M_udp_c' in s else \
                            self.model.getReaction(r).getReagent(st_ump).getCoefficient() if 'M_ump_c' in s else 1
                    c += c_utp

                if 'M_ctp_c' in s:
                    # Calculate CTP consumption/production
                    c_ctp = self.model.getReaction(r).getReagent(st_ctp).getCoefficient()
                    c_cxp = self.model.getReaction(r).getReagent(st_cdp).getCoefficient() if 'M_cdp_c' in s else \
                            self.model.getReaction(r).getReagent(st_cmp).getCoefficient() if 'M_cmp_c' in s else 1
                    c += c_ctp


                # Accumulate total consumption and production based on flux direction
                if c < 0 and active_fluxdict[r] > 0:
                    total_consumption += np.abs(c * active_fluxdict[r])
                elif c > 0 and active_fluxdict[r] > 0:
                    total_production += np.abs(c * active_fluxdict[r])
                elif c < 0 and active_fluxdict[r] < 0:
                    total_production += np.abs(c * active_fluxdict[r])
                elif c > 0 and active_fluxdict[r] < 0:
                    total_consumption += np.abs(c * active_fluxdict[r])

        return [total_consumption, total_production] 

    def generate_NADH_list(self):
        """
        Generates a list of reactions in the model that are associated with NADH.

        Returns
        -------
        reacs : list
            List of unique reaction IDs that involve NADH.
        """
        NADH_reacs = self.model.getReactionIdsAssociatedWithSpecies('M_nadh_c')
        
        # Extract reaction IDs
        reacs = [r[0] for r in NADH_reacs]
        return list(set(reacs))

    def generate_NADPH_list(self):
        """
        Generates a list of reactions in the model that are associated with NADPH.

        Returns
        -------
        reacs : list
            List of unique reaction IDs that involve NADPH.
        """
        NADPH_reacs = self.model.getReactionIdsAssociatedWithSpecies('M_nadph_c')
        
        # Extract reaction IDs
        reacs = [r[0] for r in NADPH_reacs]
        return list(set(reacs))
    
    def calculate_NADH_prod(self, active_fluxdict, split_names=True):
        """
        Calculates the total NADH production and consumption from active reactions in the flux dictionary.
    
        Parameters
        ----------
        active_fluxdict : dict
            Dictionary containing reaction IDs as keys and their corresponding flux values as values.
        split_names : bool, optional
            Whether to use split names for reactions when determining NADH and NAD+ coefficients. 
            Default is True.
    
        Returns
        -------
        list
            A list containing:
            - Total NADH production
            - Total NADH consumption
        """
        tot_nadh_prod = 0
        tot_nadh_cons = 0

        # Iterate over the list of NADH-associated reactions
        for r in self.NADH_list:
            if r in active_fluxdict:
                # Get the list of species IDs involved in the reaction
                s = self.model.getReaction(r).getSpeciesIds()
    
                # Determine the appropriate species IDs based on reaction direction and split_names flag
                if r[-4:] == '_fwd' and split_names:
                    st_nadh = r + '_M_nadh_c_fwd'
                    st_nad = r + '_M_nad_c_fwd'
                elif r[-4:] == '_rev' and split_names:
                    st_nadh = r + '_M_nadh_c_rev'
                    st_nad = r + '_M_nad_c_rev'
                else:
                    st_nadh = r + '_M_nadh_c'
                    st_nad = r + '_M_nad_c'
    
                # Initialize coefficients for NADH and NAD+
                c = 0
                c_nad = 0
                
                if 'M_nadh_c' in s:
                    # Get the coefficient for NADH and NAD+ in the reaction
                    c_nadh = self.model.getReaction(r).getReagent(st_nadh).getCoefficient()
                    
                    if 'M_nad_c' in s:
                        c_nad = self.model.getReaction(r).getReagent(st_nad).getCoefficient()
                    
                    # Determine the effective coefficient for NADH based on NADH and NAD+ availability
                    if np.abs(c_nadh) <= np.abs(c_nad): 
                        # NADH is equal to or less than NAD+, so use NADH coefficient
                        c = c_nadh
                    else:
                        # NADH is more than NAD+, so use NAD+ coefficient
                        c = c_nad
    
                # Calculate NADH production or consumption based on reaction flux and coefficient
                if c > 0 and active_fluxdict[r] > 0:
                    tot_nadh_prod += np.abs(c) * active_fluxdict[r]
                elif c < 0 and active_fluxdict[r] > 0:
                    tot_nadh_cons += np.abs(c) * active_fluxdict[r]
                elif c > 0 and active_fluxdict[r] < 0:
                    tot_nadh_cons += np.abs(c) * active_fluxdict[r]
                elif c < 0 and active_fluxdict[r] < 0:
                    tot_nadh_prod += np.abs(c) * active_fluxdict[r]
    
        return [tot_nadh_prod, tot_nadh_cons] 
    
    
    def calculate_NADPH_prod(self, active_fluxdict, split_names=True):
        """
        Calculates the total NADPH production and consumption from active reactions in the flux dictionary.
    
        Parameters
        ----------
        active_fluxdict : dict
            Dictionary containing reaction IDs as keys and their corresponding flux values as values.
        split_names : bool, optional
            Whether to use split names for reactions when determining NADPH and NADP+ coefficients. 
            Default is True.
    
        Returns
        -------
        list
            A list containing:
            - Total NADPH production
            - Total NADPH consumption
        """
        # Initialize counters for NADPH production and consumption
        tot_nadph_prod = 0
        tot_nadph_cons = 0

        # Iterate over each NADPH-associated reaction
        for r in self.NADPH_list:
            if r in active_fluxdict:
                # Get the list of species IDs involved in the reaction
                s = self.model.getReaction(r).getSpeciesIds()
                
                # Determine the appropriate species IDs based on reaction direction and split_names flag
                if r[-4:] == '_fwd' and split_names:
                    st_nadph = r + '_M_nadph_c_fwd'
                    st_nadp = r + '_M_nadp_c_fwd'
                elif r[-4:] == '_rev' and split_names:
                    st_nadph = r + '_M_nadph_c_rev'
                    st_nadp = r + '_M_nadp_c_rev'
                else:
                    st_nadph = r + '_M_nadph_c'
                    st_nadp = r + '_M_nadp_c'
    
                # Initialize coefficients for NADPH and NADP+
                c = 0
                c_nadp = 0
                
                if 'M_nadph_c' in s:
                    # Get the coefficient for NADPH and NADP+ in the reaction
                    c_nadph = self.model.getReaction(r).getReagent(st_nadph).getCoefficient()
                    
                    if 'M_nadp_c' in s:
                        c_nadp = self.model.getReaction(r).getReagent(st_nadp).getCoefficient()
                    
                    # Determine the effective coefficient for NADPH based on NADPH and NADP+ availability
                    if np.abs(c_nadph) <= np.abs(c_nadp): 
                        # If NADPH is equal to or less than NADP+, use NADPH coefficient
                        c = c_nadph
                    else:
                        # If NADPH is more than NADP+, use NADP+ coefficient
                        c = c_nadp
    
                # Calculate NADPH production or consumption based on reaction flux and coefficient
                if c > 0:
                    # Positive coefficient means production of NADPH
                    tot_nadph_prod += np.abs(c) * active_fluxdict[r]
                elif c < 0:
                    # Negative coefficient means consumption of NADPH
                    tot_nadph_cons += np.abs(c) * active_fluxdict[r]
                elif c > 0 and active_fluxdict[r] < 0:
                    # If coefficient is positive but flux is negative, this indicates NADPH consumption
                    tot_nadph_cons += np.abs(c) * active_fluxdict[r]
                elif c < 0 and active_fluxdict[r] > 0:
                    # If coefficient is negative but flux is positive, this indicates NADPH production
                    tot_nadph_prod += np.abs(c) * active_fluxdict[r]
    
        return [tot_nadph_prod, tot_nadph_cons] 