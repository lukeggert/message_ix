* Standalone call example:
* gams MACRO/macro_household_stage2_mcp.gms --in="<macro-data.gdx>" --stage1="<macro-results.gdx>" --out="<household-stage2-mcp.gdx>" --region="R12_AFR"

$TITLE MACRO Household Stage 2 MCP
$EOLCOM #
$IF NOT SET in $ABORT "no input data file provided! Use --in=<file>."
$IF NOT EXIST '%in%' $ABORT "input GDX file '%in%' does not exist!"
$IF NOT SET stage1 $ABORT "no stage1 results GDX file provided! Use --stage1=<file>."
$IF NOT EXIST '%stage1%' $ABORT "stage1 results GDX file '%stage1%' does not exist!"
$IF NOT SET out $SETGLOBAL out "output/MacroHouseholdStage2MCP.gdx"
$IF NOT SET region $SETGLOBAL region "R12_AFR"

* Household KKT system based directly on the representative-agent utility
* and the capital equation with investment substituted out.
* All households share the same functional form and exogenous prices.
* The only heterogeneity imposed here is the initial capital endowment.

SETS
    node                                 'Nodes'
    year_all                             'Full model horizon'
    sector                               'Energy sectors'
    quantile                             'Household income quantiles' / D1*D10 /
    node_macro(node)                     'Nodes active in MACRO'
    macro_horizon(year_all)              'Periods in the MACRO horizon including base period'
    seq_period(year_all, year_all)       'Mapping from one period to the next'
    last_period(year_all)                'Last model period'
    macro_base_period(year_all)          'MACRO base period'
;

ALIAS (year_all, year, year2);
ALIAS (quantile, quantile2);

SCALARS
    pareto_alpha                         'Pareto shape parameter for income distribution' / 2.0 /
    hh_scale                             'Utility scaling' / 1000 /
    pi_r                                 'Return premium by quantile rank' / 0.05 /
    kap_floor_scale                      'Minimum fraction of aggregate capital tied to target shares' / 0.70 /
    kap_floor_shape                      'Tilt of the target capital share toward upper deciles' / 1.35 /
;

PARAMETERS
    duration_period(year_all)            'Duration of one period in years'
    depr(node)                           'Annual depreciation rate'
    lotol(node)                          'Tolerance factor for lower bounds'
    EMIN(node)                           'Minimum energy requirement'
    k0(node)                             'Base-year capital'
    c0(node)
    alpha(node)                          'Aggregate non-energy utility share'
    beta_rc_spec(node)                   'Aggregate residential specific energy share'
    beta_rc_therm(node)                  'Aggregate residential thermal energy share'
    beta_transport(node)                 'Aggregate transport energy share'
    labor(node, year_all)                'Aggregate labor path'
    udf(node, year_all)                  'Utility discount factor'
    finite_time_corr(node, year_all)     'Finite horizon correction factor'
    eneprice(node, sector, year_all)     'Energy price path used in stage 2'
    quantile_share(quantile)             'Income share by quantile'
    kap_share_target(quantile)           'Target capital share profile for lower bounds'
    return_weight(quantile)              'Normalized decile-specific return weight'
    utility_weight(node)                 'Common household utility weight'
    cons_cost_factor(node)               'Consumption wedge from substituted I equation'
    labor_q_exo(node, year_all, quantile)'Exogenous labor allocation by quantile'
    energy_floor_cost_exo(node, year_all, quantile) 'Mandatory energy expenditure block'
    interest_q_exo(node, year_all, quantile) 'Decile-specific interest rate preserving the stage-1 mean'
    ret_q_exo(node, year_all, quantile)  'Effective decile-specific return term'
    ret_exo(node, year_all)              'Effective return term from stage 1 interest path'
    kgrow_exo(node, year_all)            'Exogenous terminal capital growth from stage 1'
    wage_exo(node, year_all)             'Exogenous wage path from stage 1'
    c_ref(node, year_all)                'Aggregate consumption path from stage 1'
    k_ref(node, year_all)                'Aggregate capital path from stage 1'
    c_hh_agg(node, year_all)             'Aggregate household consumption for reporting'
    k_hh_agg(node, year_all)             'Aggregate household capital for reporting'
;

POSITIVE VARIABLES
    C_REF_VAR(node, year_all)            'Aggregate consumption variable loaded from stage 1'
    K_REF_VAR(node, year_all)            'Aggregate capital variable loaded from stage 1'
    INTEREST_REF(node, year_all)         'Aggregate interest path loaded from stage 1'
    KGROW_REF(node, year_all)            'Aggregate capital growth path loaded from stage 1'
    WAGE_REF(node, year_all)             'Aggregate wage path loaded from stage 1'
    CON_HH(node, year_all, quantile)     'Non-energy household consumption by quantile'
    KAP_HH(node, year_all, quantile)     'Household capital by quantile'
;

FREE VARIABLES
    MU_CON_AGG(node, year_all)           'Shadow value of aggregate consumption matching'
    MU_KAP_AGG(node, year_all)           'Shadow value of aggregate capital matching'
    MU_CAPITAL(node, year_all, quantile) 'Capital-law multiplier'
    NU_TERM(node, year_all, quantile)    'Terminal-condition multiplier'
;

VARIABLES
    HH_UTILITY_REP                        'Reporting only: sum of household utilities'
;

EQUATIONS
    AGG_CON_MATCH(node, year_all)          'Aggregate consumption must match stage 1'
    AGG_KAP_MATCH(node, year_all)          'Aggregate capital must match stage 1'
    CAPITAL_HH(node, year_all, quantile)   'Household capital accumulation with I substituted out'
    TERMINAL_HH(node, year_all, quantile)  'Household terminal condition with I substituted out'
    FOC_CON(node, year_all, quantile)      'FOC for household consumption'
    FOC_KAP(node, year_all, quantile)      'FOC for household capital'
    HH_UTILITY_REP_DEF                     'Reporting equation only'
;

$GDXIN '%in%'
$LOAD node
$LOAD year_all = year
$LOAD sector
$LOAD node_macro
$LOAD macro_horizon
$LOAD seq_period
$LOAD last_period
$LOAD macro_base_period
$LOAD duration_period
$LOAD depr
$LOAD lotol
$LOAD EMIN
$LOAD k0,c0
$LOAD alpha
$LOAD beta_rc_spec
$LOAD beta_rc_therm
$LOAD beta_transport
$LOAD labor
$LOAD udf
$LOAD finite_time_corr
$LOAD eneprice
$GDXIN

node_macro(node)$(NOT sameas(node, '%region%')) = no ;

if (SUM(node_macro, 1) = 0,
    abort 'Selected region not found in node_macro set: %region%'
) ;

quantile_share(quantile) =
    (1 - (ORD(quantile) - 1) / CARD(quantile)) ** ((pareto_alpha - 1) / pareto_alpha)
    - (1 - ORD(quantile) / CARD(quantile)) ** ((pareto_alpha - 1) / pareto_alpha) ;

kap_share_target(quantile) = quantile_share(quantile) ** kap_floor_shape ;
kap_share_target(quantile) = kap_share_target(quantile) / SUM(quantile2, kap_share_target(quantile2)) ;

return_weight(quantile) = (1 + pi_r) ** ORD(quantile) ;
return_weight(quantile) = CARD(quantile) * return_weight(quantile) / SUM(quantile2, return_weight(quantile2)) ;

utility_weight(node_macro) = alpha(node_macro) + beta_rc_spec(node_macro) + beta_rc_therm(node_macro) + beta_transport(node_macro) ;
cons_cost_factor(node_macro) = utility_weight(node_macro) / alpha(node_macro) ;

execute_load '%stage1%',
    C_REF_VAR = C,
    K_REF_VAR = K,
    INTEREST_REF = INTEREST,
    KGROW_REF = KGROW,
    WAGE_REF = WAGE ;

c_ref(node_macro, macro_horizon) = C_REF_VAR.L(node_macro, macro_horizon) ;
k_ref(node_macro, macro_horizon) = K_REF_VAR.L(node_macro, macro_horizon) ;
ret_exo(node_macro, macro_horizon) =
    ((1 + INTEREST_REF.L(node_macro, macro_horizon)) ** duration_period(macro_horizon) - 1)
    / duration_period(macro_horizon) ;
ret_q_exo(node_macro, macro_horizon, quantile) =
    ret_exo(node_macro, macro_horizon) * return_weight(quantile) ;
interest_q_exo(node_macro, macro_horizon, quantile) =
    (1 + ret_q_exo(node_macro, macro_horizon, quantile) * duration_period(macro_horizon)) ** (1 / duration_period(macro_horizon)) - 1 ;
kgrow_exo(node_macro, macro_horizon) = KGROW_REF.L(node_macro, macro_horizon) ;
wage_exo(node_macro, macro_horizon) = WAGE_REF.L(node_macro, macro_horizon) ;

* Household labor income is heterogeneous and follows quantile income shares.
labor_q_exo(node_macro, macro_horizon, quantile) = labor(node_macro, macro_horizon) * quantile_share(quantile) ;

energy_floor_cost_exo(node_macro, macro_horizon, quantile) =
    labor_q_exo(node_macro, macro_horizon, quantile) * EMIN(node_macro) * (
          eneprice(node_macro, 'rc_spec', macro_horizon) / 1000
        + eneprice(node_macro, 'rc_therm', macro_horizon) / 1000
        + eneprice(node_macro, 'transport', macro_horizon) / 1000
    ) ;

CON_HH.LO(node_macro, macro_horizon, quantile)$(NOT macro_base_period(macro_horizon)) = lotol(node_macro) * c0(node_macro) * quantile_share(quantile) ;
KAP_HH.LO(node_macro, macro_base_period, quantile) = lotol(node_macro) * k0(node_macro) * quantile_share(quantile) ;
KAP_HH.LO(node_macro, macro_horizon, quantile)$(NOT macro_base_period(macro_horizon)) =
    MAX(
        lotol(node_macro) * k0(node_macro) * quantile_share(quantile),
        kap_floor_scale * k_ref(node_macro, macro_horizon) * kap_share_target(quantile)
    ) ;

*CON_HH.L(node_macro, macro_horizon, quantile)$(NOT macro_base_period(macro_horizon)) = lotol(node_macro) * quantile_share(quantile) ;
KAP_HH.L(node_macro, macro_base_period, quantile) = k0(node_macro) * quantile_share(quantile) ;
KAP_HH.L(node_macro, macro_horizon, quantile)$(NOT macro_base_period(macro_horizon)) =
    MAX(
        k0(node_macro) * quantile_share(quantile),
        kap_floor_scale * k_ref(node_macro, macro_horizon) * kap_share_target(quantile)
    ) ;

KAP_HH.FX(node_macro, macro_base_period, quantile) = k0(node_macro) * quantile_share(quantile) ;
CON_HH.FX(node_macro, macro_base_period, quantile) = 0 ;

AGG_CON_MATCH(node_macro, year)$(macro_horizon(year) AND NOT macro_base_period(year))..
    c_ref(node_macro, year) =E= SUM(quantile, CON_HH(node_macro, year, quantile)) ;

AGG_KAP_MATCH(node_macro, year)$(macro_horizon(year) AND NOT macro_base_period(year))..
    k_ref(node_macro, year) =E= SUM(quantile, KAP_HH(node_macro, year, quantile)) ;

CAPITAL_HH(node_macro, year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year))..
KAP_HH(node_macro, year, quantile) =E=
    SUM(year2$seq_period(year2, year),
        KAP_HH(node_macro, year2, quantile)
            * ((1 - depr(node_macro)) ** duration_period(year) + duration_period(year) * ret_q_exo(node_macro, year2, quantile))
        + duration_period(year) * (
            labor_q_exo(node_macro, year2, quantile) * wage_exo(node_macro, year2)
            - energy_floor_cost_exo(node_macro, year2, quantile)
            - cons_cost_factor(node_macro) * CON_HH(node_macro, year2, quantile)
        )
    ) ;

TERMINAL_HH(node_macro, last_period, quantile)..
    KAP_HH(node_macro, last_period, quantile)
        * (ret_q_exo(node_macro, last_period, quantile) - kgrow_exo(node_macro, last_period) - depr(node_macro))
    + labor_q_exo(node_macro, last_period, quantile) * wage_exo(node_macro, last_period)
    - energy_floor_cost_exo(node_macro, last_period, quantile)
    - cons_cost_factor(node_macro) * CON_HH(node_macro, last_period, quantile)
=E= 0 ;

FOC_CON(node_macro, year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year))..
    hh_scale * udf(node_macro, year) * utility_weight(node_macro) * duration_period(year)
        / CON_HH(node_macro, year, quantile)
    - MU_CON_AGG(node_macro, year)
    - SUM(year2$seq_period(year, year2),
        MU_CAPITAL(node_macro, year2, quantile) * duration_period(year2) * cons_cost_factor(node_macro)
    )
    - NU_TERM(node_macro, year, quantile)$last_period(year) * cons_cost_factor(node_macro)
    =G= 0 ;

FOC_KAP(node_macro, year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year))..
    MU_CAPITAL(node_macro, year, quantile)
    - MU_KAP_AGG(node_macro, year)
    - SUM(year2$seq_period(year, year2),
        MU_CAPITAL(node_macro, year2, quantile)
            * ((1 - depr(node_macro)) ** duration_period(year2) + duration_period(year2) * ret_q_exo(node_macro, year, quantile))
    )
    - NU_TERM(node_macro, year, quantile)$last_period(year)
        * (ret_q_exo(node_macro, year, quantile) - kgrow_exo(node_macro, year) - depr(node_macro))
    =G= 0 ;

HH_UTILITY_REP_DEF..
HH_UTILITY_REP =E=
SUM(node_macro,
    hh_scale * (
        SUM((year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year) AND NOT last_period(year)),
            udf(node_macro, year) * (
                utility_weight(node_macro) * LOG(CON_HH(node_macro, year, quantile))
                - beta_rc_spec(node_macro) * LOG(eneprice(node_macro, 'rc_spec', year) / 1000)
                + beta_rc_spec(node_macro) * LOG(beta_rc_spec(node_macro) / alpha(node_macro))
                - beta_rc_therm(node_macro) * LOG(eneprice(node_macro, 'rc_therm', year) / 1000)
                + beta_rc_therm(node_macro) * LOG(beta_rc_therm(node_macro) / alpha(node_macro))
                - beta_transport(node_macro) * LOG(eneprice(node_macro, 'transport', year) / 1000)
                + beta_transport(node_macro) * LOG(beta_transport(node_macro) / alpha(node_macro))
            ) * duration_period(year)
        )
        + SUM((year, quantile)$last_period(year),
            udf(node_macro, year) * (
                utility_weight(node_macro) * LOG(CON_HH(node_macro, year, quantile))
                - beta_rc_spec(node_macro) * LOG(eneprice(node_macro, 'rc_spec', year) / 1000)
                + beta_rc_spec(node_macro) * LOG(beta_rc_spec(node_macro) / alpha(node_macro))
                - beta_rc_therm(node_macro) * LOG(eneprice(node_macro, 'rc_therm', year) / 1000)
                + beta_rc_therm(node_macro) * LOG(beta_rc_therm(node_macro) / alpha(node_macro))
                - beta_transport(node_macro) * LOG(eneprice(node_macro, 'transport', year) / 1000)
                + beta_transport(node_macro) * LOG(beta_transport(node_macro) / alpha(node_macro))
            ) * duration_period(year)
            + 1 / finite_time_corr(node_macro, year)
        )
    )
) ;

MODEL MACRO_HOUSEHOLD_STAGE2_MCP /
    AGG_CON_MATCH.MU_CON_AGG
    AGG_KAP_MATCH.MU_KAP_AGG
    CAPITAL_HH.MU_CAPITAL
    TERMINAL_HH.NU_TERM
    FOC_CON.CON_HH
    FOC_KAP.KAP_HH
/ ;

OPTION MCP = PATH ;

SOLVE MACRO_HOUSEHOLD_STAGE2_MCP USING MCP ;

c_hh_agg(node_macro, macro_horizon) = SUM(quantile, CON_HH.L(node_macro, macro_horizon, quantile)) ;
k_hh_agg(node_macro, macro_horizon) = SUM(quantile, KAP_HH.L(node_macro, macro_horizon, quantile)) ;

execute_unload '%out%',
    CON_HH,
    KAP_HH,
    MU_CON_AGG,
    MU_KAP_AGG,
    MU_CAPITAL,
    NU_TERM,
    c_ref,
    k_ref,
    c_hh_agg,
    k_hh_agg,
    interest_q_exo,
    ret_q_exo,
    ret_exo,
    wage_exo,
    kgrow_exo,
    labor_q_exo,
    energy_floor_cost_exo,
    quantile_share,
    HH_UTILITY_REP ;