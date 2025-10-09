POSITIVE VARIABLES
    K(node, year_all)                Capital stock in period year
    Y(node, year_all)                Production in period year

    PHYSENE(node, sector, year_all)  Physical end-use service or commodity use
    PRODENE(node, sector, year_all)  Value of end-use services or commodities in the production function
    NEWENE(node, sector, year_all)   New end-use service or commodity (production function value)

    C(node, year_all)                Consumption (Trillion $)
    I(node, year_all)                Investment (Trillion $)
;

VARIABLES
    UTILITY                          Utility function (discounted log of consumption)
    EC(node, year_all)               System costs (Trillion $) based on MESSAGE model run
;

EQUATIONS
    UTILITY_FUNCTION                      Utility function (discounted log of consumption summed over all projection periods)
    CAPITAL_CONSTRAINT(node, year_all)    Capital constraint

    CAPITAL(node, year_all)           New capital

    INVESTMENT(node, year_all)            Total capital stock across all vintages
    PRODUCTION(node, year_all)      Total production across all vintages

    NEW_ENERGY(node, sector, year_all)    New end-use services or commodities (production function)
    ENERGY_SUPPLY(node, sector, *)        Supply of end-use services or commodities

    COST_ENERGY(node, year_all)           system costs approximation based on MESSAGE input
    TERMINAL_CONDITION(node, year_all)    Terminal condition
;

UTILITY_FUNCTION..
UTILITY =E=
SUM(node_active,
    1000 * (SUM(year $ (NOT macro_base_period(year) AND NOT last_period(year)),
            udf(node_active, year) * ( (alpha(node_active) + beta1(node_active)) * LOG(C(node_active, year)) - beta1(node_active) * LOG(eneprice(node_active, 'light', year)) + beta1(node_active) * LOG(beta1(node_active)/alpha(node_active)) ) * duration_period(year) )
        + SUM(year $ last_period(year),
            udf(node_active, year) * ( (alpha(node_active) + beta1(node_active)) * LOG(C(node_active, year)) - beta1(node_active) * LOG(eneprice(node_active, 'light', year)) + beta1(node_active) * LOG(beta1(node_active)/alpha(node_active)) ) * (duration_period(year) ) + 1/finite_time_corr(node_active, year)) )
)
;

CAPITAL_CONSTRAINT(node_active, year)..
Y(node_active, year) =E=
C(node_active, year) + I(node_active, year) + EC(node_active, year)
;

CAPITAL(node_active, year) $ (NOT macro_base_period(year))..
K(node_active, year) =E=
SUM(year2$( seq_period(year2,year) ), K(node_active, year2) * (1 - depr(node_active))**duration_period(year) + I(node_active, year2))
;

INVESTMENT(node_active, year) $ (NOT macro_base_period(year))..
I(node_active, year) =E= K(node_active, year) * ((1 + interestrate(year))**duration_period(year) - 1) + duration_period(year) * labor(node_active, year) * wage(node_active) - duration_period(year) * eneprice(node_active, 'light', year) * labor(node_active, year) * EMIN(node_active) - duration_period(year) * (((alpha(node_active) + beta1(node_active))/alpha(node_active))) * C(node_active, year)
;

PRODUCTION(node_active, year) $ (NOT macro_base_period(year))..
Y(node_active, year) =E=
(LAKL(node_active) * K(node_active, year)**(rho(node_active) * kpvs(node_active)) * newlab(node_active, year)**(RHO(node_active) * (1 - kpvs(node_active)))
+ PRFCONST(node_active, 'heat') * PRODENE(node_active, 'heat', year)**rho(node_active)) **(1/rho(node_active))
;

NEW_ENERGY(node_active, sector, year) $ (NOT macro_base_period(year))..
PRODENE(node_active, sector, year) =E=
SUM(year2$( seq_period(year2,year) ), PRODENE(node_active, sector, year2)) * (1 - depr(node_active))**duration_period(year) + NEWENE(node_active, sector, year)
;

ENERGY_SUPPLY(node_active, sector, year) $ (NOT macro_base_period(year))..
PHYSENE(node_active, sector, year) =G=
PRODENE(node_active, sector, year) * aeei_factor(node_active, sector, year)
;

COST_ENERGY(node_active, year) $ (NOT macro_base_period(year))..
EC(node_active, year) =E=
(total_cost(node_active, year)/1000
+ SUM(sector, eneprice(node_active, sector, year) * 1E-3 * (PHYSENE(node_active, sector, year) - enestart(node_active, sector, year)))
+ SUM(sector, eneprice(node_active, sector, year) * 1E-3 / enestart(node_active, sector, year) * (PHYSENE(node_active, sector, year) - enestart(node_active, sector, year)) * (PHYSENE(node_active, sector, year) - enestart(node_active, sector, year))))
;

TERMINAL_CONDITION(node_active, last_period)..
I(node_active, last_period) =G= K(node_active, last_period) * (grow(node_active, last_period) + depr(node_active))
;

* ------------------------------------------------------------------------------
* model definition
* ------------------------------------------------------------------------------

MODEL MESSAGE_MACRO /
    UTILITY_FUNCTION
    CAPITAL_CONSTRAINT
    CAPITAL
    PRODUCTION
    INVESTMENT
    NEW_ENERGY
    ENERGY_SUPPLY
    COST_ENERGY
    TERMINAL_CONDITION
/ ;

MESSAGE_MACRO.optfile = 1;
