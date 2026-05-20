* Standalone call example:
* gams MACRO/macro_household_stage2.gms --in="<macro-data.gdx>" --stage1="<macro-results.gdx>" --out="<household-stage2.gdx>"

$TITLE MACRO Household Stage 2
$EOLCOM #
$IF NOT SET in $ABORT "no input data file provided! Use --in=<file>."
$IF NOT EXIST '%in%' $ABORT "input GDX file '%in%' does not exist!"
$IF NOT SET stage1 $ABORT "no stage1 results GDX file provided! Use --stage1=<file>."
$IF NOT EXIST '%stage1%' $ABORT "stage1 results GDX file '%stage1%' does not exist!"
$IF NOT SET out $SETGLOBAL out "output/MacroHouseholdStage2.gdx"
$IF NOT SET region $SETGLOBAL region "R12_AFR"

* Household Stage 2 (decentralized quantile block with exogenous macro paths)
*
* Purpose:
* - Run after the aggregate MACRO model has been solved.
* - Take prices and aggregate labor paths as given.
* - Solve independent forward-looking household problems for each quantile.
* - Report implied aggregate C and K from household decisions.
*
* Important:
* - This file is fully standalone: no include files are required.
* - All required sets, parameters and stage-1 variables are loaded here from GDX.
* - Aggregate consistency with stage 1 is checked by reporting only; it is not
*   enforced as a hard constraint here.

SETS
    node                                 'Nodes'
    type_node                            'Node categories'
    cat_node(type_node, node)            'Mapping of node categories to nodes'
    year_all                             'Full model horizon'
    type_year                            'Year categories'
    cat_year(type_year, year_all)        'Mapping of years to categories'
    sector                               'Energy sectors'
    quantile                             'Household income quantiles' / D1*D10 /
    node_macro(node)                     'Nodes active in MACRO'
    historical(year_all)                 'Periods before the model horizon'
    model_horizon(year_all)              'Periods in the model horizon'
    macro_horizon(year_all)              'Periods in the MACRO horizon including base period'
    seq_period(year_all, year_all)       'Mapping from one period to the next'
    map_period(year_all, year_all)       'Mapping from one period to itself and all successors'
    first_period(year_all)               'First model period'
    last_period(year_all)                'Last model period'
    macro_initial_period(year_all)       'Initialization period for MACRO'
    macro_base_period(year_all)          'MACRO base period'
;

ALIAS (year_all, year, year2, year_all2, year_all3);

SCALARS
    pareto_alpha                         'Pareto shape parameter for income distribution' / 2.0 /
    hh_scale                             'Objective scaling' / 1000 /
;

PARAMETERS
    duration_period(year_all)            'Duration of one period in years'
    depr(node)                           'Annual depreciation rate'
    lotol(node)                          'Tolerance factor for lower bounds'
    EMIN(node)                           'Minimum energy requirement'
    k_final(node)                        'Terminal aggregate capital target'
    alpha(node)                          'Aggregate non-energy utility share'
    beta_rc_spec(node)                   'Aggregate residential specific energy share'
    beta_rc_therm(node)                  'Aggregate residential thermal energy share'
    beta_transport(node)                 'Aggregate transport energy share'

    labor(node, year_all)                'Aggregate labor path'
    udf(node, year_all)                  'Utility discount factor'
    finite_time_corr(node, year_all)     'Finite-horizon correction term'
    eneprice(node, sector, year_all)     'Energy price path used in stage 2'
    gdp_base(node)                       'Base-year GDP'
    i0(node)                             'Base-year investment'
    c0(node)                             'Base-year consumption'
    k0(node)                             'Base-year capital'
    quantile_share(quantile)             'Income share by quantile'
    utility_weight(node)                 'Common household utility weight'
    cons_cost_factor(node)               'Common household consumption wedge'

    interest_exo(node, year_all)         'Exogenous interest path from stage 1'
    wage_exo(node, year_all)             'Exogenous wage path from stage 1'
    kgrow_exo(node, year_all)            'Exogenous capital growth path from stage 1'
    labor_exo(node, year_all)            'Exogenous aggregate labor path'
    labor_q_exo(node, year_all, quantile)'Exogenous labor allocation by quantile'
    energy_floor_cost_exo(node, year_all, quantile) 'Mandatory energy expenditure block'
    k_macro_ref(node, year_all)          'Aggregate K path from stage 1 for reporting only'
    c_macro_ref(node, year_all)          'Aggregate C path from stage 1 for reporting only'
    kap_init(node, quantile)             'Initial capital by quantile in base period'
    kap_terminal_target(node, quantile)  'Terminal capital target by quantile'
    utility_weight_q(node, quantile)     'Utility weight for each quantile'
    cons_cost_factor_q(node, quantile)   'Consumption cost wedge in budget law'
    c_hh_gap(node, year_all)             'Household aggregate C minus stage-1 C'
    k_hh_gap(node, year_all)             'Household aggregate K minus stage-1 K'
;

POSITIVE VARIABLES
    K(node, year_all)                    'Aggregate capital path loaded from stage 1'
    C(node, year_all)                    'Aggregate consumption path loaded from stage 1'
    I(node, year_all)                    'Aggregate investment path loaded from stage 1'
    CON_HH(node, year_all, quantile)     'Non-energy household consumption by quantile'
    I_HH(node, year_all, quantile)       'Household investment by quantile'
    KAP_HH(node, year_all, quantile)     'Household capital by quantile'
    C_HH_AGG(node, year_all)             'Aggregate household consumption implied by stage 2'
    K_HH_AGG(node, year_all)             'Aggregate household capital implied by stage 2'
;

VARIABLES
    INTEREST(node, year_all)             'Interest path loaded from stage 1'
    WAGE(node, year_all)                 'Wage path loaded from stage 1'
    KGROW(node, year_all)                'Capital growth path loaded from stage 1'
    HH_UTILITY                           'Sum of separable household utilities'
;

EQUATIONS
    HH_UTILITY_FUNCTION                  'Stage-2 household utility'
    HH_CAPITAL_EVOLUTION(node, year_all, quantile) 'Household capital accumulation'
    HH_INVESTMENT(node, year_all, quantile) 'Household investment accounting'
    HH_TERMINAL_CONDITION(node, year_all, quantile) 'Household terminal condition'
    HH_AGG_CON_DEF(node, year_all)       'Aggregate household consumption definition'
    HH_AGG_KAP_DEF(node, year_all)       'Aggregate household capital definition'
;

* ------------------------------------------------------------------------------
* Load all required symbols directly from the input GDX.
* ------------------------------------------------------------------------------

$GDXIN '%in%'
$LOAD node
$LOAD type_node
$LOAD cat_node
$LOAD year_all = year
$LOAD type_year
$LOAD cat_year
$LOAD sector
$LOAD node_macro
$LOAD historical
$LOAD model_horizon
$LOAD macro_horizon
$LOAD seq_period
$LOAD map_period
$LOAD first_period
$LOAD last_period
$LOAD macro_initial_period
$LOAD macro_base_period
$LOAD duration_period
$LOAD depr
$LOAD lotol
$LOAD EMIN
$LOAD k_final
$LOAD alpha
$LOAD beta_rc_spec
$LOAD beta_rc_therm
$LOAD beta_transport
$LOAD labor
$LOAD udf
$LOAD finite_time_corr
$LOAD eneprice
$LOAD gdp_base
$LOAD i0
$LOAD c0
$LOAD k0
$GDXIN

node_macro(node)$(NOT sameas(node, '%region%')) = no ;

if (SUM(node_macro, 1) = 0,
    abort 'Selected region not found in node_macro set: %region%'
) ;

* ------------------------------------------------------------------------------
* The only object rebuilt locally is the quantile income share.
* ------------------------------------------------------------------------------

quantile_share(quantile) =
    (1 - (ORD(quantile) - 1) / CARD(quantile)) ** ((pareto_alpha - 1) / pareto_alpha)
    - (1 - ORD(quantile) / CARD(quantile)) ** ((pareto_alpha - 1) / pareto_alpha) ;

utility_weight(node_macro) = alpha(node_macro) + beta_rc_spec(node_macro) + beta_rc_therm(node_macro) + beta_transport(node_macro) ;
cons_cost_factor(node_macro) = utility_weight(node_macro) / alpha(node_macro) ;

* ------------------------------------------------------------------------------
* Load the solved aggregate paths from the stage-1 results GDX.
* ------------------------------------------------------------------------------

execute_load '%stage1%', K, C, I, INTEREST, WAGE, KGROW ;

interest_exo(node_macro, macro_horizon) = INTEREST.L(node_macro, macro_horizon) ;
wage_exo(node_macro, macro_horizon) = WAGE.L(node_macro, macro_horizon) ;
kgrow_exo(node_macro, macro_horizon) = KGROW.L(node_macro, macro_horizon) ;
labor_exo(node_macro, macro_horizon) = labor(node_macro, macro_horizon) ;
labor_q_exo(node_macro, macro_horizon, quantile) = labor_exo(node_macro, macro_horizon) * quantile_share(quantile) ;

energy_floor_cost_exo(node_macro, macro_horizon, quantile) =
    labor_q_exo(node_macro, macro_horizon, quantile) * EMIN(node_macro) * (
          eneprice(node_macro, 'rc_spec', macro_horizon) / 1000
        + eneprice(node_macro, 'rc_therm', macro_horizon) / 1000
        + eneprice(node_macro, 'transport', macro_horizon) / 1000
    ) ;

k_macro_ref(node_macro, macro_horizon) = K.L(node_macro, macro_horizon) ;
c_macro_ref(node_macro, macro_horizon) = C.L(node_macro, macro_horizon) ;

kap_init(node_macro, quantile) = k0(node_macro) * quantile_share(quantile) ;
kap_terminal_target(node_macro, quantile) = k_final(node_macro) * quantile_share(quantile) ;

utility_weight_q(node_macro, quantile) = utility_weight(node_macro) ;
cons_cost_factor_q(node_macro, quantile) = cons_cost_factor(node_macro) ;

* ------------------------------------------------------------------------------
* Bounds and initial values.
* ------------------------------------------------------------------------------

CON_HH.LO(node_macro, macro_horizon, quantile) = lotol(node_macro) * c0(node_macro) * quantile_share(quantile) ;
I_HH.LO(node_macro, macro_horizon, quantile) = 0 ;
KAP_HH.LO(node_macro, macro_horizon, quantile) = lotol(node_macro) * k0(node_macro) * quantile_share(quantile) ;

CON_HH.L(node_macro, macro_horizon, quantile) = C.L(node_macro, macro_horizon) * quantile_share(quantile) ;
I_HH.L(node_macro, macro_horizon, quantile) = I.L(node_macro, macro_horizon) * quantile_share(quantile) ;
KAP_HH.L(node_macro, macro_horizon, quantile) = K.L(node_macro, macro_horizon) * quantile_share(quantile) ;

KAP_HH.FX(node_macro, macro_base_period, quantile) = kap_init(node_macro, quantile) ;
KAP_HH.FX(node_macro, last_period, quantile) = kap_terminal_target(node_macro, quantile) ;

* ------------------------------------------------------------------------------
* Objective and constraints.
* ------------------------------------------------------------------------------

HH_UTILITY_FUNCTION..
HH_UTILITY =E=
SUM(node_macro,
    hh_scale * (
        SUM(year$(macro_horizon(year) AND NOT macro_base_period(year) AND NOT last_period(year)),
            udf(node_macro, year)
            * duration_period(year)
            * SUM(quantile,
                utility_weight(node_macro) * LOG(CON_HH(node_macro, year, quantile))
            )
        )
        + SUM(year$last_period(year),
            udf(node_macro, year) * (
                duration_period(year)
                * SUM(quantile,
                    utility_weight(node_macro) * LOG(CON_HH(node_macro, year, quantile))
                )
                + 1 / finite_time_corr(node_macro, year)
            )
        )
    )
) ;

HH_CAPITAL_EVOLUTION(node_macro, year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year))..
KAP_HH(node_macro, year, quantile) =E=
    SUM(year2$seq_period(year2, year),
        KAP_HH(node_macro, year2, quantile) * (1 - depr(node_macro)) ** duration_period(year)
        + duration_period(year) * I_HH(node_macro, year2, quantile)
    ) ;

HH_INVESTMENT(node_macro, year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year))..
I_HH(node_macro, year, quantile) =E=
    (KAP_HH(node_macro, year, quantile)
        * ((1 + interest_exo(node_macro, year)) ** duration_period(year) - 1)
        / duration_period(year))
    + labor_q_exo(node_macro, year, quantile) * wage_exo(node_macro, year)
    - energy_floor_cost_exo(node_macro, year, quantile)
    - cons_cost_factor(node_macro) * CON_HH(node_macro, year, quantile) ;

HH_TERMINAL_CONDITION(node_macro, last_period, quantile)..
I_HH(node_macro, last_period, quantile) =E=
    KAP_HH(node_macro, last_period, quantile) * (kgrow_exo(node_macro, last_period) + depr(node_macro)) ;

HH_AGG_CON_DEF(node_macro, year)$macro_horizon(year)..
C_HH_AGG(node_macro, year) =E= SUM(quantile, CON_HH(node_macro, year, quantile)) ;

HH_AGG_KAP_DEF(node_macro, year)$macro_horizon(year)..
K_HH_AGG(node_macro, year) =E= SUM(quantile, KAP_HH(node_macro, year, quantile)) ;

MODEL MACRO_HOUSEHOLD_STAGE2 /
    HH_UTILITY_FUNCTION
    HH_CAPITAL_EVOLUTION
    HH_INVESTMENT
    HH_TERMINAL_CONDITION
    HH_AGG_CON_DEF
    HH_AGG_KAP_DEF
/ ;

MACRO_HOUSEHOLD_STAGE2.optfile = 0 ;

SOLVE MACRO_HOUSEHOLD_STAGE2 USING NLP MAXIMIZING HH_UTILITY ;

abort$(MACRO_HOUSEHOLD_STAGE2.modelstat = 4 OR MACRO_HOUSEHOLD_STAGE2.modelstat = 5)
    'Household stage 2 did not solve to a bounded optimum.' ;

c_hh_gap(node_macro, macro_horizon) = C_HH_AGG.L(node_macro, macro_horizon) - c_macro_ref(node_macro, macro_horizon) ;
k_hh_gap(node_macro, macro_horizon) = K_HH_AGG.L(node_macro, macro_horizon) - k_macro_ref(node_macro, macro_horizon) ;

DISPLAY HH_UTILITY.L, C_HH_AGG.L, K_HH_AGG.L, c_hh_gap, k_hh_gap, CON_HH.L, KAP_HH.L ;

execute_unload '%out%',
    CON_HH,
    I_HH,
    KAP_HH,
    C_HH_AGG,
    K_HH_AGG,
    c_hh_gap,
    k_hh_gap,
    interest_exo,
    wage_exo,
    kgrow_exo,
    labor_exo,
    labor_q_exo,
    energy_floor_cost_exo,
    k_macro_ref,
    c_macro_ref,
    kap_init,
    kap_terminal_target,
    quantile_share,
    utility_weight,
    cons_cost_factor,
    HH_UTILITY ;
