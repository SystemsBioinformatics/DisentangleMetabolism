import numpy as np
import pandas as pd
import re
import cbmpy

class PrecursorSeparator:
    def __init__(self, model, result_df, biomass_reac='R_EX_BIOMASS'):
        self.model = model
        self.result_df = result_df
        self.biomass_reac = biomass_reac

        # Internal state
        self.catabolic_intermediates = []
        self.catabolic_external = []
        self.energy_carriers = []
        self.added_supplies = []
        self.ex_added = []
        self.kos = []
        self.fluxes = None
        self.mceq = None
        self.fva_df = None

        # Define possible energy carriers once
        self.possible_energy_carriers = [
            'M_atp_c', 'M_adp_c', 'M_nad_c', 'M_nadh_c', 'M_nadph_c', 'M_nadp_c',
            'M_fdxr_42_c', 'M_fdxo_42_c', 'M_fad_c', 'M_fadh2_c', 'M_q8_c', 'M_q8h2_c',
            'M_atp_m', 'M_adp_m', 'M_nad_m', 'M_nadh_m', 'M_nadph_m', 'M_nadp_m',
            'M_f420_DASH_2h2_c','M_f420_DASH_2_c', 'M_fdred_c', 'M_fdox_c',
            'M_mphenh2_c','M_mphen_c','M_omchr_e','M_omcho_e','M_focytC_c','M_ficytC_c',
            'M_flxso_c','M_flxr_c','M_mql8_c','M_mqn8_c','M_focytc_m','M_ficytc_m',
            'M_q6h2_m','M_q6_m','M_fad_m','M_fadh2_m','M_lpam_m','M_dhlam_m',
            'M_gdp_c','M_gtp_c','M_ttp_c','M_tdp_c','M_utp_c','M_udp_c',
            'M_fdxo_2_2_c','M_fdxr_c','M_nmn_c','M_nmnh_c']

    def process_inputs(self):
        catabolic_intermediates = []
        catabolic_external = []
        energy_carriers = []

        for r in self.result_df.index:
            if np.abs(self.result_df.loc[r, 'catabolic']) > 1e-5:
                if 'demand' in r:
                    continue

                species = self.model.getReaction(r).getSpeciesIds()
                for s in species:
                    if s[-1] == 'e' and s not in catabolic_external and s not in self.possible_energy_carriers:
                        catabolic_external.append(s)
                for s in species:
                    if s not in catabolic_intermediates and (s[-1] == 'c' or s[-1] == 'm') and s not in self.possible_energy_carriers:
                        catabolic_intermediates.append(s)
                    elif s not in energy_carriers and s in self.possible_energy_carriers:
                        energy_carriers.append(s)

        if 'M_h_e' not in catabolic_external:
            catabolic_external.append('M_h_e')
        if 'M_h_c' not in catabolic_intermediates:
            catabolic_intermediates.append('M_h_c')

        self.catabolic_intermediates = catabolic_intermediates
        self.catabolic_external = catabolic_external
        self.energy_carriers = energy_carriers
        return catabolic_intermediates, catabolic_external, energy_carriers

    def create_species_and_reactions(self):
        ex_added = []
        for s in self.catabolic_intermediates:
            if s in self.energy_carriers:
                continue

            name = self.model.getSpecies(s).name
            idd = s[:-1] + 'e'
            if self.model.getSpecies(idd) is None:
                self.model.createSpecies(idd, name=name + ' ex', compartment='e')

        for s in self.catabolic_intermediates:
            self.create_transport_reaction(s)
            r_name = self.create_exchange_reaction(s)
            if r_name:
                ex_added.append(r_name)

        self.ex_added = ex_added
        return ex_added

    def create_transport_reaction(self, s):
        self.model.createReaction(f'R_{s[2:]}_trans', name=f'{s[2:-2]} transport', reversible=True)
        self.model.createReactionReagent(f'R_{s[2:]}_trans', s, -1)
        self.model.createReactionReagent(f'R_{s[2:]}_trans', s[:-1] + 'e', 1)
        self.model.setReactionBounds(f'R_{s[2:]}_trans', -1000, 1000)

    def create_exchange_reaction(self, s):
        name = self.model.getSpecies(s).name
        chemform = self.model.getSpecies(s).chemFormula
        match = re.search(r'P(\d*)', chemform)
        number_after_P = int(match.group(1)) if match and match.group(1) else 1

        r_ex_st = 'R_EX_' + s[2:-2] + '_e'
        if r_ex_st not in self.model.getReactionIds() and r_ex_st + '_rev' not in self.model.getReactionIds() and r_ex_st+'_fwd' not in self.model.getReactionIds():
            self.model.createReaction(r_ex_st, name=s[2:-2] + ' exchange', reversible=True)
            self.model.createReactionReagent(r_ex_st, s[:-1]+'e', -1)
            if 'CoA' in name and s != 'M_coa_c':
                self.model.createReactionReagent(r_ex_st, 'M_coa_e', 1)
            elif 'thf' in s and s != 'M_thf_c':
                if 'M_thf_e' not in self.model.getSpeciesIds():
                    self.model.createSpecies('M_thf_e')
                self.model.createReactionReagent(r_ex_st, 'M_thf_e', 1)
            elif ('phosph' in name or 'Phosph' in name) and s != 'M_pi_c':
                self.model.createReactionReagent(r_ex_st, 'M_pi_e', number_after_P)
            self.model.setReactionBounds(r_ex_st, -1000, 1000)
            return r_ex_st
        return None

    def add_supply_reactions(self):
        # Supply reactions defined with metabolite IDs and coefficients
        supplies = {
            'NADPH': [('M_nadph_c', 1.0), ('M_nadp_c', -1.0), ('M_h_c', -1.0)],
            'NADH': [('M_nadh_c', 1.0), ('M_nad_c', -1.0), ('M_h_c', -1.0)],
            'ATP': [('M_atp_c', 1.0), ('M_h2o_c', 1.0), ('M_adp_c', -1.0), ('M_h_c', -1.0), ('M_pi_c', -1.0)],
            'FDXR': [('M_fdxo_42_c', -1.0), ('M_fdxr_42_c', 1.0)],
            'FADH': [('M_fadh2_c', 1.0), ('M_fad_c', -1.0), ('M_h_c', -2.0)],
            'Q8H2': [('M_q8h2_c', 1.0), ('M_q8_c', -1.0)],
            'ATP_m': [('M_atp_m', 1.0), ('M_h2o_m', 1.0), ('M_adp_m', -1.0), ('M_h_m', -1.0), ('M_pi_m', -1.0)],
            'NADH_m': [('M_nadh_m', 1.0), ('M_nad_m', -1.0), ('M_h_m', -1.0)],
            'NADPH_m': [('M_nadph_m', 1.0), ('M_nadp_m', -1.0), ('M_h_m', -1.0)],
            'F420R': [('M_f420_DASH_2h2_c', 1.0), ('M_f420_DASH_2_c', -1.0)],
            'FDRED':[('M_fdred_c', 1.0), ('M_fdox_c', -1.0)],
            'FDXR2': [('M_fdxrd_c', 1.0), ('M_fdxo_2_2_c', -1.0)],
            'MPHEN': [('M_mphenh2_c', 1.0), ('M_mphen_c', -1.0)],
            'OMCHR': [('M_omchr_e', 1.0), ('M_omcho_e', -1.0)],
            'MQN8': [('M_mqn8_c', 1.0), ('M_mql8_c', -1.0)],
            'FLXR': [('M_flxr_c', 1.0), ('M_flxso_c', -1.0)],
            'FICYTC': [('M_ficytC_c', 1.0), ('M_focytC_c', -1.0)],
            'FICYTC_m':[('M_focytc_m', -1.0), ('M_ficytc_m', 1.0)],
            'Q6H2_m':[('M_q6h2_m', 1.0), ('M_q6_m', -1.0)],
            'FADH2_m':[('M_fadh2_m', 1.0), ('M_fad_m', -1.0)],
            'DHLAM_m':[('M_dhlam_m', 1.0), ('M_lpam_m', -1.0)],
            'GTP': [('M_gtp_c', 1.0), ('M_h2o_c', 1.0), ('M_gdp_c', -1.0), ('M_h_c', -1.0), ('M_pi_c', -1.0)],
            'UTP': [('M_utp_c', 1.0), ('M_h2o_c', 1.0), ('M_udp_c', -1.0), ('M_h_c', -1.0), ('M_pi_c', -1.0)]
        }

        added_supplies = []
        for name, reagents in supplies.items():
            if any(metabolite in self.energy_carriers for metabolite, _ in reagents):
                self.model.createReaction(f'R_{name}_supply', reversible=True, name=f'{name} supply')
                for metabolite, coefficient in reagents:
                    self.model.createReactionReagent(f'R_{name}_supply', metabolite=metabolite, coefficient=coefficient)
                self.model.setReactionBounds(f'R_{name}_supply', -1000, 1000)
                added_supplies.append(f'R_{name}_supply')

        self.added_supplies = added_supplies
        return added_supplies

    def adjust_reaction_bounds(self):
        kos = []
        for m in self.catabolic_external:
            rid = 'R_EX' + m[1:]
            if rid in self.model.getReactionIds():
                self.model.setReactionBounds(rid, -1000, 1000)
                self.model.getReaction(rid).reversible = True
            elif rid+'_fwd' in self.model.getReactionIds():
                self.model.setReactionBounds(rid+'_fwd', -1000, 1000)
            elif rid+'_rev' in self.model.getReactionIds():
                self.model.setReactionBounds(rid+'_rev', -1000, 1000)
            else:
                raise Exception('Missing exchange reaction?')

        for r in self.result_df.index:
            if 'demand' in r or 'ATPS' in r:
                continue
            catabolic_reac = False
            if np.abs(self.result_df.loc[r, 'catabolic']) > 1e-5:
                catabolic_reac = True
                species = self.model.getReaction(r).getSpeciesIds()
                for s in species:
                    if s not in self.catabolic_intermediates and s not in self.energy_carriers:
                        catabolic_reac = False
            if catabolic_reac and 'ADK1' not in r:
                self.model.setReactionBounds(r, 0, 0)
                kos.append(r)
        self.kos = kos
        return kos

    def get_MCEQ_precursors(self, obj, fluxes):
        mceq_model = self.model.clone()
        fluxes_mceq = fluxes.copy()
        for r in self.added_supplies:
            print(r)
            mceq_model.deleteReactionAndBounds(r+'_rev')
            mceq_model.deleteReactionAndBounds(r+'_fwd')
            fluxes_mceq.pop(r+'_rev')
            fluxes_mceq.pop(r+'_fwd')
        mceq = obj.get_MCEQ(mceq_model, biomass_reac=self.biomass_reac,
                            flux_dist=pd.Series(fluxes_mceq),
                            additional_active_metabolites=self.energy_carriers)
        return mceq

    def run(self, obj):
        self.process_inputs()
        self.create_species_and_reactions()
        self.add_supply_reactions()
        self.adjust_reaction_bounds()

        self.model.setObjectiveFlux(self.biomass_reac, osense='minimize')
        model_notsplit = self.model.clone()
        
        self.model = cbmpy.CBTools.splitReversibleReactions(self.model)

        result = cbmpy.CBGLPK.glpk_analyzeModel(self.model, method='s')
        self.fluxes = self.model.getReactionValues()
        self.mceq = self.get_MCEQ_precursors(obj, self.fluxes)

        
        result = cbmpy.CBGLPK.glpk_analyzeModel(model_notsplit, method='s')
        fva = cbmpy.CBCPLEX.cplx_FluxVariabilityAnalysis(model_notsplit, selected_reactions=self.ex_added+self.added_supplies)

        self.fva_df = pd.DataFrame(fva[0], index=fva[1],
                                   columns=['Reaction', 'Reduced Costs', 'Variability Min', 'Variability Max', 'abs(Max-Min)', 'MinStatus', 'MaxStatus'])

        return self.fluxes, self.mceq, self.fva_df, self.kos
