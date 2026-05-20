* Standalone call example:
* gams /home/lukas/environments/macro_uba/lib/python3.10/site-packages/message_ix/model/MACRO/macro_household_clean_mcp_yfixed.gms --in="/home/lukas/environments/macro_uba/lib/python3.10/site-packages/message_ix/model/output/MsgOutput_MESSAGEix_ssp2_baseline_2304_add_macro_price_10_loob_eneprice.gdx" --stage1="/home/lukas/environments/macro_uba/lib/python3.10/site-packages/message_ix/model/output/MsgOutput_MESSAGEix_ssp2_baseline_2304_add_macro_price_10_loob_eneprice.gdx" --out="/home/lukas/household_clean_mcp_price.gdx" --region="R12_WEU"

$TITLE MACRO Household Clean MCP With Fixed Y Path
$EOLCOM #
$IF NOT SET in $ABORT "no input data file provided! Use --in=<file>."
$IF NOT EXIST '%in%' $ABORT "input GDX file '%in%' does not exist!"
$IF NOT SET stage1 $ABORT "no stage1 results GDX file provided! Use --stage1=<file>."
$IF NOT EXIST '%stage1%' $ABORT "stage1 results GDX file '%stage1%' does not exist!"
$IF NOT SET out $SETGLOBAL out "output/MacroHouseholdCleanMCP.gdx"
$IF NOT SET region $SETGLOBAL region "R12_AFR"

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

SCALARS
    pareto_alpha                         'Pareto shape parameter for income distribution' / 2.0 /
    hh_scale                             'Utility scaling' / 1000 /
;

PARAMETERS
    duration_period(year_all)            'Duration of one period in years'
    depr(node)                           'Annual depreciation rate'
    lotol(node)                          'Tolerance factor for lower bounds'
    EMIN(node)                           'Minimum energy requirement'
    k0(node)                             'Base-year capital'
    c0(node)                             'Base-year consumption'
    alpha(node)                          'Aggregate non-energy utility share'
    beta_rc_spec(node)                   'Aggregate residential specific energy share'
    beta_rc_therm(node)                  'Aggregate residential thermal energy share'
    beta_transport(node)                 'Aggregate transport energy share'
    labor(node, year_all)                'Aggregate labor path'
    udf(node, year_all)                  'Utility discount factor'
    finite_time_corr(node, year_all)     'Finite horizon correction factor'
    eneprice(node, sector, year_all)     'Energy price path'
    kpvs(node)                           'Capital value share parameter'
    esub(node)                           'Elasticity of substitution'
    lakl(node)                           'Production function coefficient of capital and labor'
    rho(node)                            'Production function exponent'
    quantile_share(quantile)             'Income share by quantile'
    utility_weight(node)                 'Common household utility weight'
    cons_cost_factor(node)               'Consumption wedge from substituted I equation'
    labor_q_exo(node, year_all, quantile)'Exogenous labor allocation by quantile'
    energy_floor_cost_exo(node, year_all, quantile) 'Mandatory energy expenditure block'
    y_ref(node, year_all)                'Stage-1 output path used as fixed firm output'
    c_ref(node, year_all)                'Stage-1 aggregate consumption path'
    k_ref(node, year_all)                'Stage-1 aggregate capital path'
    wage_ref(node, year_all)             'Stage-1 wage path'
    interest_ref(node, year_all)         'Stage-1 interest path'
    kgrow_ref(node, year_all)            'Stage-1 terminal capital growth path'
    E(quantile, node, sector, year_all)   'Energy consumption by quantile, sector and year (post-solve reporting)'
    c_mcp_agg(node, year_all)            'Aggregate MCP consumption for reporting'
    k_mcp_agg(node, year_all)            'Aggregate MCP capital for reporting'
    k_gap(node, year_all)                'Difference MCP minus stage-1 capital'
    c_gap(node, year_all)                'Difference MCP minus stage-1 consumption'
    wage_gap(node, year_all)             'Difference MCP minus stage-1 wage'
    interest_gap(node, year_all)         'Difference MCP minus stage-1 interest'
    k_gap_rel(node, year_all)            'Relative capital difference'
    c_gap_rel(node, year_all)            'Relative consumption difference'
    wage_gap_rel(node, year_all)         'Relative wage difference'
    interest_gap_rel(node, year_all)     'Relative interest difference'
;

POSITIVE VARIABLES
    Y_REF_VAR(node, year_all)            'Stage-1 output variable loaded from GDX'
    C_REF_VAR(node, year_all)            'Stage-1 consumption variable loaded from GDX'
    K_REF_VAR(node, year_all)            'Stage-1 capital variable loaded from GDX'
    INTEREST_REF_VAR(node, year_all)     'Stage-1 interest variable loaded from GDX'
    WAGE_REF_VAR(node, year_all)         'Stage-1 wage variable loaded from GDX'
    KGROW_REF_VAR(node, year_all)        'Stage-1 capital growth variable loaded from GDX'
    CON_HH(node, year_all, quantile)     'Non-energy household consumption by quantile'
    KAP_HH(node, year_all, quantile)     'Household capital by quantile'
    K_FIRM(node, year_all)               'Firm capital demand'
    L_FIRM(node, year_all)               'Firm labor demand'
    INTEREST_EQ(node, year_all)          'Equilibrium interest rate'
    WAGE_EQ(node, year_all)              'Equilibrium wage rate'
;

FREE VARIABLES
    MU_CAPITAL(node, year_all, quantile) 'Household capital-law multiplier'
    NU_TERM(node, year_all, quantile)    'Household terminal-condition multiplier'
;

VARIABLES
    HH_UTILITY_REP                        'Reporting only: sum of household utilities'
;

EQUATIONS
    CAPITAL_MARKET(node, year_all)        'Capital market clearing'
    LABOR_MARKET(node, year_all)          'Labor market clearing'
    FIRM_FOC_CAPITAL(node, year_all)      'Firm capital demand condition with fixed Y path'
    FIRM_FOC_LABOR(node, year_all)        'Firm labor demand condition with fixed Y path'
    CAPITAL_HH(node, year_all, quantile)  'Household capital accumulation with price-taking factor returns'
    TERMINAL_HH(node, year_all, quantile) 'Household terminal condition'
    FOC_CON(node, year_all, quantile)     'FOC for household consumption'
    FOC_KAP(node, year_all, quantile)     'FOC for household capital'
    HH_UTILITY_REP_DEF                    'Reporting equation only'
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
$LOAD kpvs,esub,lakl
$GDXIN

node_macro(node)$(NOT sameas(node, '%region%')) = no ;

if (SUM(node_macro, 1) = 0,
    abort 'Selected region not found in node_macro set: %region%'
) ;

rho(node_macro) = (esub(node_macro) - 1) / esub(node_macro) ;

quantile_share(quantile) =
    (1 - (ORD(quantile) - 1) / CARD(quantile)) ** ((pareto_alpha - 1) / pareto_alpha)
    - (1 - ORD(quantile) / CARD(quantile)) ** ((pareto_alpha - 1) / pareto_alpha) ;

utility_weight(node_macro) = alpha(node_macro) + beta_rc_spec(node_macro) + beta_rc_therm(node_macro) + beta_transport(node_macro) ;
cons_cost_factor(node_macro) = utility_weight(node_macro) / alpha(node_macro) ;

execute_load '%stage1%',
    Y_REF_VAR = Y,
    C_REF_VAR = C,
    K_REF_VAR = K,
    INTEREST_REF_VAR = INTEREST,
    WAGE_REF_VAR = WAGE,
    KGROW_REF_VAR = KGROW ;

y_ref(node_macro, macro_horizon) = Y_REF_VAR.L(node_macro, macro_horizon) ;
c_ref(node_macro, macro_horizon) = C_REF_VAR.L(node_macro, macro_horizon) ;
k_ref(node_macro, macro_horizon) = K_REF_VAR.L(node_macro, macro_horizon) ;
interest_ref(node_macro, macro_horizon) = INTEREST_REF_VAR.L(node_macro, macro_horizon) ;
wage_ref(node_macro, macro_horizon) = WAGE_REF_VAR.L(node_macro, macro_horizon) ;
kgrow_ref(node_macro, macro_horizon) = KGROW_REF_VAR.L(node_macro, macro_horizon) ;

labor_q_exo(node_macro, macro_horizon, quantile) = labor(node_macro, macro_horizon) * quantile_share(quantile) ;

energy_floor_cost_exo(node_macro, macro_horizon, quantile) =
    labor_q_exo(node_macro, macro_horizon, quantile) * EMIN(node_macro) * (
          eneprice(node_macro, 'rc_spec', macro_horizon) / 1000
        + eneprice(node_macro, 'rc_therm', macro_horizon) / 1000
        + eneprice(node_macro, 'transport', macro_horizon) / 1000
    ) ;

CON_HH.LO(node_macro, macro_horizon, quantile)$(NOT macro_base_period(macro_horizon)) = lotol(node_macro) * c0(node_macro) * quantile_share(quantile) ;
KAP_HH.LO(node_macro, macro_horizon, quantile) = lotol(node_macro) * k0(node_macro) * quantile_share(quantile) ;
K_FIRM.LO(node_macro, macro_horizon)$(NOT macro_base_period(macro_horizon)) = lotol(node_macro) * k0(node_macro) ;
L_FIRM.LO(node_macro, macro_horizon)$(NOT macro_base_period(macro_horizon)) = lotol(node_macro) * labor(node_macro, macro_horizon) ;
INTEREST_EQ.LO(node_macro, macro_horizon)$(NOT macro_base_period(macro_horizon)) = lotol(node_macro) ;
WAGE_EQ.LO(node_macro, macro_horizon)$(NOT macro_base_period(macro_horizon)) = lotol(node_macro) ;

CON_HH.L(node_macro, macro_horizon, quantile)$(NOT macro_base_period(macro_horizon)) = MAX(lotol(node_macro) * c0(node_macro) * quantile_share(quantile), c_ref(node_macro, macro_horizon) * quantile_share(quantile)) ;
KAP_HH.L(node_macro, macro_horizon, quantile) = MAX(lotol(node_macro) * k0(node_macro) * quantile_share(quantile), k_ref(node_macro, macro_horizon) * quantile_share(quantile)) ;
K_FIRM.L(node_macro, macro_horizon) = MAX(lotol(node_macro) * k0(node_macro), k_ref(node_macro, macro_horizon)) ;
L_FIRM.L(node_macro, macro_horizon) = MAX(lotol(node_macro) * labor(node_macro, macro_horizon), labor(node_macro, macro_horizon)) ;
INTEREST_EQ.L(node_macro, macro_horizon) = MAX(lotol(node_macro), interest_ref(node_macro, macro_horizon)) ;
WAGE_EQ.L(node_macro, macro_horizon) = MAX(lotol(node_macro), wage_ref(node_macro, macro_horizon)) ;

KAP_HH.FX(node_macro, macro_base_period, quantile) = k0(node_macro) * quantile_share(quantile) ;
CON_HH.FX(node_macro, macro_base_period, quantile) = 0 ;
K_FIRM.FX(node_macro, macro_base_period) = k_ref(node_macro, macro_base_period) ;
L_FIRM.FX(node_macro, macro_base_period) = labor(node_macro, macro_base_period) ;
INTEREST_EQ.FX(node_macro, macro_base_period) = MAX(lotol(node_macro), interest_ref(node_macro, macro_base_period)) ;
WAGE_EQ.FX(node_macro, macro_base_period) = MAX(lotol(node_macro), wage_ref(node_macro, macro_base_period)) ;

CAPITAL_MARKET(node_macro, year)$(macro_horizon(year) AND NOT macro_base_period(year))..
    SUM(quantile, KAP_HH(node_macro, year, quantile)) - K_FIRM(node_macro, year) =G= 0 ;

LABOR_MARKET(node_macro, year)$(macro_horizon(year) AND NOT macro_base_period(year))..
    labor(node_macro, year) - L_FIRM(node_macro, year) =G= 0 ;

FIRM_FOC_CAPITAL(node_macro, year)$(macro_horizon(year) AND NOT macro_base_period(year))..
    INTEREST_EQ(node_macro, year) =G=
        y_ref(node_macro, year) ** (1 - rho(node_macro))
        * lakl(node_macro) * kpvs(node_macro)
        * K_FIRM(node_macro, year) ** (rho(node_macro) * kpvs(node_macro) - 1)
        * L_FIRM(node_macro, year) ** (rho(node_macro) * (1 - kpvs(node_macro))) ;

FIRM_FOC_LABOR(node_macro, year)$(macro_horizon(year) AND NOT macro_base_period(year))..
    WAGE_EQ(node_macro, year) =G=
        y_ref(node_macro, year) ** (1 - rho(node_macro))
        * lakl(node_macro) * (1 - kpvs(node_macro))
        * K_FIRM(node_macro, year) ** (rho(node_macro) * kpvs(node_macro))
        * L_FIRM(node_macro, year) ** (rho(node_macro) * (1 - kpvs(node_macro)) - 1) ;

CAPITAL_HH(node_macro, year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year))..
    KAP_HH(node_macro, year, quantile) =E=
        SUM(year2$seq_period(year2, year),
            KAP_HH(node_macro, year2, quantile) * (1 - depr(node_macro)) ** duration_period(year)
            + duration_period(year) * (
                KAP_HH(node_macro, year2, quantile)
                    * (((1 + INTEREST_EQ(node_macro, year2)) ** duration_period(year2) - 1) / duration_period(year2))
                + labor_q_exo(node_macro, year2, quantile) * WAGE_EQ(node_macro, year2)
                - energy_floor_cost_exo(node_macro, year2, quantile)
                - cons_cost_factor(node_macro) * CON_HH(node_macro, year2, quantile)
            )
        ) ;

TERMINAL_HH(node_macro, last_period, quantile)..
    KAP_HH(node_macro, last_period, quantile)
        * ((((1 + INTEREST_EQ(node_macro, last_period)) ** duration_period(last_period) - 1) / duration_period(last_period))
            - kgrow_ref(node_macro, last_period) - depr(node_macro))
    + labor_q_exo(node_macro, last_period, quantile) * WAGE_EQ(node_macro, last_period)
    - energy_floor_cost_exo(node_macro, last_period, quantile)
    - cons_cost_factor(node_macro) * CON_HH(node_macro, last_period, quantile)
=E= 0 ;

FOC_CON(node_macro, year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year))..
    hh_scale * udf(node_macro, year) * utility_weight(node_macro) * duration_period(year)
        / CON_HH(node_macro, year, quantile)
    - SUM(year2$seq_period(year, year2),
        MU_CAPITAL(node_macro, year2, quantile) * duration_period(year2) * cons_cost_factor(node_macro)
    )
    - NU_TERM(node_macro, year, quantile)$last_period(year) * cons_cost_factor(node_macro)
    =G= 0 ;

FOC_KAP(node_macro, year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year))..
    MU_CAPITAL(node_macro, year, quantile)
    - SUM(year2$seq_period(year, year2),
        MU_CAPITAL(node_macro, year2, quantile)
            * ((1 - depr(node_macro)) ** duration_period(year2)
            + duration_period(year2) * (((1 + INTEREST_EQ(node_macro, year)) ** duration_period(year2) - 1) / duration_period(year2)))
    )
    - NU_TERM(node_macro, year, quantile)$last_period(year)
        * ((((1 + INTEREST_EQ(node_macro, year)) ** duration_period(year) - 1) / duration_period(year))
            - kgrow_ref(node_macro, year) - depr(node_macro))
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

MODEL MACRO_HOUSEHOLD_CLEAN_MCP_YFIXED /
    CAPITAL_MARKET.INTEREST_EQ
    LABOR_MARKET.WAGE_EQ
    FIRM_FOC_CAPITAL.K_FIRM
    FIRM_FOC_LABOR.L_FIRM
    CAPITAL_HH.MU_CAPITAL
    TERMINAL_HH.NU_TERM
    FOC_CON.CON_HH
    FOC_KAP.KAP_HH
/ ;

OPTION MCP = PATH ;

SOLVE MACRO_HOUSEHOLD_CLEAN_MCP_YFIXED USING MCP ;

* Energy consumption by quantile - ENERGY_ACCOUNTING2
E(quantile, node_macro, 'rc_spec', year)$(macro_horizon(year) AND NOT macro_base_period(year)) =
    EMIN(node_macro)
    + (beta_rc_spec(node_macro) / alpha(node_macro))
    * CON_HH.L(node_macro, year, quantile)
    / (eneprice(node_macro, 'rc_spec', year) / 1000) ;

E(quantile, node_macro, 'rc_therm', year)$(macro_horizon(year) AND NOT macro_base_period(year)) =
    EMIN(node_macro)
    + (beta_rc_therm(node_macro) / alpha(node_macro))
    * CON_HH.L(node_macro, year, quantile)
    / (eneprice(node_macro, 'rc_therm', year) / 1000) ;

E(quantile, node_macro, 'transport', year)$(macro_horizon(year) AND NOT macro_base_period(year)) =
    EMIN(node_macro)
    + (beta_transport(node_macro) / alpha(node_macro))
    * CON_HH.L(node_macro, year, quantile)
    / (eneprice(node_macro, 'transport', year) / 1000) ;

c_mcp_agg(node_macro, macro_horizon) = SUM(quantile, CON_HH.L(node_macro, macro_horizon, quantile)) ;
k_mcp_agg(node_macro, macro_horizon) = SUM(quantile, KAP_HH.L(node_macro, macro_horizon, quantile)) ;

k_gap(node_macro, macro_horizon) = k_mcp_agg(node_macro, macro_horizon) - k_ref(node_macro, macro_horizon) ;
c_gap(node_macro, macro_horizon) = c_mcp_agg(node_macro, macro_horizon) - c_ref(node_macro, macro_horizon) ;
wage_gap(node_macro, macro_horizon) = WAGE_EQ.L(node_macro, macro_horizon) - wage_ref(node_macro, macro_horizon) ;
interest_gap(node_macro, macro_horizon) = INTEREST_EQ.L(node_macro, macro_horizon) - interest_ref(node_macro, macro_horizon) ;

k_gap_rel(node_macro, macro_horizon)$k_ref(node_macro, macro_horizon) = 100 * k_gap(node_macro, macro_horizon) / k_ref(node_macro, macro_horizon) ;
c_gap_rel(node_macro, macro_horizon)$c_ref(node_macro, macro_horizon) = 100 * c_gap(node_macro, macro_horizon) / c_ref(node_macro, macro_horizon) ;
wage_gap_rel(node_macro, macro_horizon)$wage_ref(node_macro, macro_horizon) = 100 * wage_gap(node_macro, macro_horizon) / wage_ref(node_macro, macro_horizon) ;
interest_gap_rel(node_macro, macro_horizon)$interest_ref(node_macro, macro_horizon) = 100 * interest_gap(node_macro, macro_horizon) / interest_ref(node_macro, macro_horizon) ;

execute_unload '%out%',
    CON_HH,
    KAP_HH,
    K_FIRM,
    L_FIRM,
    INTEREST_EQ,
    WAGE_EQ,
    MU_CAPITAL,
    NU_TERM,
    y_ref,
    c_ref,
    k_ref,
    wage_ref,
    interest_ref,
    kgrow_ref,
    c_mcp_agg,
    k_mcp_agg,
    k_gap,
    c_gap,
    wage_gap,
    interest_gap,
    k_gap_rel,
    c_gap_rel,
    wage_gap_rel,
    interest_gap_rel,
    labor_q_exo,
    energy_floor_cost_exo,
    quantile_share,
    E,
    HH_UTILITY_REP ;