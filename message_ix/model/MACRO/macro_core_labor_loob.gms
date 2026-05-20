***
* .. _macro-core:
*
* MACRO core formulation
* ======================
*
* MACRO is a macroeconomic model maximizing the intertemporal utility function of a single representative producer-consumer
* in each node (or macro-economic region). The optimization result is a sequence of optimal savings, investment, and consumption decisions.
* The main variables of the model are the capital stock, available labor, and commodity inputs, which together determine the
* total output of an economy according to a nested constant elasticity of substitution (CES) production function. End-use service
* demands in the (commercial) demand categories of MESSAGE is determined within the model, and is consistent with commodity
* supply curves, which are inputs to the model.
*
***

*----------------------------------------------------------------------------------------------------------------------*
* Notation declaration                                                                                                 *
*----------------------------------------------------------------------------------------------------------------------*

***
* Notation declaration
* ~~~~~~~~~~~~~~~~~~~~
* The following short notation is used in the mathematical description relative to the GAMS code:
*
* ============== =============================== ===================================================================
* Math Notation  GAMS set & index notation       Description
* ============== =============================== ===================================================================
* :math:`n`      node (or node_active in loops)  spatial node corresponding to the macro-economic MESSAGE regions
* :math:`y`      year                            year (2005, 2010, 2020, ..., 2100)
* :math:`s`      sector                          sector corresponding to the (commercial) end-use demands of MESSAGE
* ============== =============================== ===================================================================
*
* A listing of all parameters used in MACRO together with a decription can be found in the table below.
*
* ================================== ================================================================================================================================
* Parameter                          Description
* ================================== ================================================================================================================================
* :math:`\text{period}_t`            Number of years in time period :math:`t` (forward diff)
* :math:`\text{total_cost}_{n,t}`    Total system costs in region :math:`n` and period :math:`t` from MESSAGE model run
* :math:`\text{enestart}_{n,s,t}`    Consumption level of (commercial) end-use services :math:`s` in region :math:`n` and period :math:`t` from MESSAGE model run
* :math:`\text{eneprice}_{n,s,t}`    Shadow prices of (commercial) end-use services :math:`s` in region :math:`n` and period :math:`t` from MESSAGE model run
* :math:`\text{E}_{min,n,s,t}`       Subsistence level of direct energy consumption (end-use service) in region :math:`n`, sector :math:`s` and period :math:`t`
* :math:`\text{h}_{n,s,t}`           Share of the direct energy consumption of the total energy production in region :math:`n`, sector :math:`s` and period :math:`t`

* :math:`\epsilon_n`                 Elasticity of substitution between capital-labor and total energy in region :math:`n`
* :math:`\rho_n`                     :math:`\epsilon - 1 / \epsilon` where :math:`\epsilon` is the elasticity of substitution in region :math:`n`
* :math:`\beta_n`                    Consumption value share parameter in region :math:`n`
* :math:`\sigma_{n,s}`               Direct energy consumption value share parameter in region :math:`n` and of sector :math:`s`
* :math:`\delta_n`                   Annual depreciation rate in region :math:`n`
* :math:`\alpha_n`                   Capital value share parameter in region :math:`n`
* :math:`a_n`                        Production function coefficient of capital and labor in region :math:`n`
* :math:`b_{n,s}`                    Production function coefficients of the different end-use sectors in region :math:`n`, sector :math:`s` and period :math:`t`
* :math:`\text{udf}_{n,t}`           Utility discount factor in period year in region :math:`n` and period :math:`t`
* :math:`\text{L}_{n,t}`             Labor force in region :math:`n` and period :math:`t`
* :math:`\text{grow}_{n,t}`          Annual growth rates of potential GDP in region :math:`n` and period :math:`t`
* :math:`\text{aeei}_{n,s,t}`        Autonomous energy efficiency improvement (AEEI) in region :math:`n`, sector :math:`s` and period :math:`t`
* :math:`\text{fin_time}_{n,t}`      Finite time horizon correction factor in utility function in region :math:`n` and period :math:`t`
* ================================== ================================================================================================================================
***

*----------------------------------------------------------------------------------------------------------------------*
* Variable definitions                                                                                                 *
*----------------------------------------------------------------------------------------------------------------------*

***
* Decision variables
* ~~~~~~~~~~~~~~~~~~~~
*
* =============================== =========================================================== ==============================================================================================================
* Variable                        Definition                                                  Description
* =============================== =========================================================== ==============================================================================================================
* :math:`\text{K}_{n,y}`          :math:`\text{K}_{n, y}\geq 0 ~ \forall n, y`                Capital stock in region :math:`n` and period :math:`y`
* :math:`\text{Y}_{n,y}`          :math:`\text{Y}_{n, y}\geq 0 ~ \forall n, y`                Total production in region :math:`n` and period :math:`y`
* :math:`\text{C}_{n,y}`          :math:`\text{C}_{n, y}\geq 0 ~ \forall n, y`                Consumption in region :math:`n` and period :math:`y`
* :math:`\text{I}_{n,y}`          :math:`\text{I}_{n, y}\geq 0 ~ \forall n, y`                Investment in region :math:`n` and period :math:`y`
* :math:`\text{WAGE}`             :math:`\text{WAGE} \in \left[-\infty..\infty\right]`        Effective wage in region :math:`n` and period :math:`y`
* :math:`\text{KGROW}`            :math:`\text{KGROW} \in \left[-\infty..\infty\right]`       Capital growth in region :math:`n` and period :math:`y`
* :math:`\text{PHYSENE}_{n,s,y}`  :math:`\text{PHYSENE}_{n, s, y}\geq 0 ~ \forall n, s, y`    Physical end-use service use in region :math:`n`, sector :math:`s` and period :math:`y`
* :math:`\text{TE}_{n,s,y}`       :math:`\text{TE}_{n, s, y}\geq 0 ~ \forall n, s, y`         Value of total end-use service in the production function and utility function in region :math:`n`, sector :math:`s` and period :math:`y`
* :math:`\text{E}_{n,s,y}`        :math:`\text{E}_{n, s, y}\geq 0 ~ \forall n, s, y`          Value of direct energy consumption of end-use service of households in the utility function in region :math:`n`, sector :math:`s` and period :math:`y`
* :math:`\text{YE}_{n,s,y}`       :math:`\text{YE}_{n, s, y}\geq 0 ~ \forall n, s, y`         Value of end-use service energy consumption in the production function in region :math:`n`, sector :math:`s` and period :math:`y`
* :math:`\text{EC}_{n,y}`         :math:`\text{EC} \in \left[-\infty..\infty\right]`          Approximation of system costs based on MESSAGE results
* :math:`\text{UTILITY}`          :math:`\text{UTILITY} \in \left[-\infty..\infty\right]`     Utility function (discounted log of consumption)
* =============================== =========================================================== ==============================================================================================================
*
***

* ------------------------------------------------------------------------------
* model variable declaration
* ------------------------------------------------------------------------------

POSITIVE VARIABLES
    K(node, year_all)                Capital stock in period year
    Y(node, year_all)                Production in period year

    PHYSENE(node, sector, year_all)  Physical end-use service or commodity use
    YE(node, sector, year_all)       Value of end-use services or commodities in the production function
    E(node, sector, year_all)        Value of end-use services or commodities in the utility function
    TE(node, sector, year_all)       Total value of end-use services or commodities

    C(node, year_all)                Consumption (Trillion $)
    I(node, year_all)                Investment (Trillion $)

    WAGE(node, year_all)             Effective wage rate 

    INTEREST(node, year_all)         Interest rate#
*AE(node, year_all)
;

VARIABLES
    UTILITY                          Utility function (discounted log of consumption)
    EC(node, year_all)               System costs (Trillion $) based on MESSAGE model run
    KGROW(node, year_all)            Capital growth
;

Variables
* auxiliary variables for demand, prices, costs and GDP (for reporting when MESSAGE is run with MACRO)
    GDP(node,year_all)               gross domestic product (GDP) in market exchange rates for MACRO reporting
;

* ------------------------------------------------------------------------------
* model equations declaration
* ------------------------------------------------------------------------------

EQUATIONS
    UTILITY_FUNCTION                      Utility function (discounted log of consumption summed over all projection periods)

    CAPITAL_CONSTRAINT(node, year_all)    Capital constraint
    CAPITAL(node, year_all)               Capital accumulation
    INVESTMENT(node, year_all)            Investment accounting

    PRODUCTION(node, year_all)            Production
*TECHNOLOGY(node, year_all) 

*TE_EQUATION(node, sector, year_all)
    ENERGY_ACCOUNTING(node, sector, year_all)  Energy accounting
    ENERGY_ACCOUNTING2(node, sector, year_all) Energy accounting (household share)
    ENERGY_SUPPLY(node, sector, *)             Supply of end-use services or commodities

    COST_ENERGY(node, year_all)           System costs approximation based on MESSAGE input
    TERMINAL_CONDITION(node, year_all)    Terminal condition

    MARGINAL_PRODUCT_CAPITAL(node, year_all)  Marginal product of capital equals interest rate plus depreciation
    MARGINAL_PRODUCT_LABOR(node, year_all)    Marginal product of labor equals wage rate
    CAPITAL_GROWTH(node, year_all)            Capital growth rate definition
;

* ------------------------------------------------------------------------------
* model equations definition
* ------------------------------------------------------------------------------

***
* Equation UTILITY_FUNCTION
* ---------------------------------
* The utility function, which is maximized, sums the discounted logarithm of consumption :math:`\text{C}_{n,y}` and direct energy consumption 
* of end-use services :math:`\text{E}_{n,s,y}` of a single representative household over the entire time horizon of the model.
*
* The utility function and the capital formulation of the optimization problem are derived in previously (See equations 10 and 13 of the model documentation). 
*
* .. math:: \text{UTILITY} = \sum_{n} \bigg( &  \sum_{y |  (  (  {ord}( y )   >  1 )  \wedge  (  {ord}( y )   <   | y |  )  )} \text{udf}_{n, y} \cdot \bigg( (\beta_n + \sum_{s=1}^{3} \sigma_{s,n}) \log(\text{C}_{n, y}) - \sum_{s=1}^{3} \sigma_{s,n} \log(\text{p}_{n,s,y}) + \sum_{s=1}^{3} \sigma_{s,n} \log\left(\frac{\sigma_{s,n}}{\beta_n}\right) \bigg) \cdot \text{duration_period}_{y} \\
*                                 + &\sum_{y |  (  {ord}( y ) =  | y | ) } \text{udf}_{n, y} \cdot \bigg( (\beta_n + \sum_{s=1}^{3} \sigma_{s,n}) \log(\text{C}_{n, y}) - \sum_{s=1}^{3} \sigma_{s,n} \log(\text{p}_{n,s,y}) + \sum_{s=1}^{3} \sigma_{s,n} \log\left(\frac{\sigma_{s,n}}{\beta_n}\right) \bigg) \cdot \big( \text{duration_period}_{y-1} + \frac{1}{\text{fin_time}_{n, y}} \big) \bigg)
*
*
***

UTILITY_FUNCTION..
UTILITY =E=
SUM(node_active,
    1000 * (SUM(year $ (NOT macro_base_period(year) AND NOT last_period(year)),
            udf(node_active, year) * ( (alpha(node_active) + beta_rc_spec(node_active) + beta_rc_therm(node_active) + beta_transport(node_active)) * LOG(C(node_active, year)) 
            - beta_rc_spec(node_active) * LOG(eneprice(node_active, 'rc_spec', year)/1000) + beta_rc_spec(node_active) * LOG(beta_rc_spec(node_active)/alpha(node_active)) 
            - beta_rc_therm(node_active) * LOG(eneprice(node_active, 'rc_therm', year)/1000) + beta_rc_therm(node_active) * LOG(beta_rc_therm(node_active)/alpha(node_active))
            - beta_transport(node_active) * LOG(eneprice(node_active, 'transport', year)/1000) + beta_transport(node_active) * LOG(beta_transport(node_active)/alpha(node_active)) ) * duration_period(year) )
        + SUM(year $ last_period(year),
            udf(node_active, year) * ( (alpha(node_active) + beta_rc_spec(node_active) + beta_rc_therm(node_active) + beta_transport(node_active)) * LOG(C(node_active, year)) 
            - beta_rc_spec(node_active) * LOG(eneprice(node_active, 'rc_spec', year)/1000) + beta_rc_spec(node_active) * LOG(beta_rc_spec(node_active)/alpha(node_active)) 
            - beta_rc_therm(node_active) * LOG(eneprice(node_active, 'rc_therm', year)/1000) + beta_rc_therm(node_active) * LOG(beta_rc_therm(node_active)/alpha(node_active))
            - beta_transport(node_active) * LOG(eneprice(node_active, 'transport', year)/1000) + beta_transport(node_active) * LOG(beta_transport(node_active)/alpha(node_active)) ) * (duration_period(year) ) + 1/finite_time_corr(node_active, year)) )
)
;

***
* Equation CAPITAL_CONSTRAINT
* ---------------------------------
* The following equation specifies the allocation of total production among current consumption :math:`\text{C}_{n, y}`, investment into building up capital stock excluding
* the sectors represented in MESSAGE :math:`\text{I}_{n, y}` and the MESSAGE system costs :math:`\text{EC}_{n, y}` which are derived from a previous MESSAGE model run. As described in :cite:`Manne-Richels-1992`, the first-order
* optimality conditions lead to the Ramsey rule for the optimal allocation of savings, investment and consumption over time.
*
* .. math:: \text{Y}_{n, y} = \text{C}_{n, y} + \text{I}_{r, y} + \text{EC}_{n, y} \qquad \forall{n, y}
*
***

CAPITAL_CONSTRAINT(node_active, year)..
Y(node_active, year) =E=
C(node_active, year) + I(node_active, year) + EC(node_active, year)
;

***
* Equation CAPITAL
* ---------------------------------
*
* The household maximizes its utility subject to the constraint on wealth accumulation of capital in the sectors not represented in the energy model MESSAGE.
* The net capital (or wealth) formation :math:`\text{K}_{n,t}` is derived from the existing capital stock, returns on capital, labor income, minus the expenses
* for direct energy consumption and all other consumption goods, as well as depreciation of the previous capital stock.
*
* .. math:: \text{K}_{t+1,n} = (1 - \delta_{n})^{\text{period}_{t}} \text{K}_{t,n} + ((1 + r_{t,n})^{\text{period}_{t}} - 1) \text{K}_{t,n} + \text{period}_{t} \left( w_{t,n} L_{t,n} - \sum_{s=1}^{3} p_{t,s,n} L_{n,t} E_{min,t,s,n} - \frac{\beta + \sum_{s=1}^{3} \sigma_{s,n}}{\beta} C_{t,n} \right) \qquad (15)
*
***

CAPITAL(node_active, year) $ (NOT macro_base_period(year))..
K(node_active, year) =E=
SUM(year2$( seq_period(year2,year) ), K(node_active, year2) * (1 - depr(node_active))**duration_period(year) + duration_period(year) * I(node_active, year2))
; 

INVESTMENT(node_active, year) $ (NOT macro_base_period(year))..
I(node_active, year) =E=
    (K(node_active, year) * ((1 + INTEREST(node_active, year))**duration_period(year) - 1) / duration_period(year))
    + labor(node_active, year) * WAGE(node_active, year) 
    - eneprice(node_active, 'rc_spec', year)/1000 * labor(node_active, year) * EMIN(node_active)
    - eneprice(node_active, 'rc_therm', year)/1000 * labor(node_active, year) * EMIN(node_active) 
    - eneprice(node_active, 'transport', year)/1000 * labor(node_active, year) * EMIN(node_active) 
    - ((alpha(node_active) + beta_rc_spec(node_active) + beta_rc_therm(node_active) + beta_transport(node_active))/alpha(node_active)) * C(node_active, year)
;

***
* Equation PRODUCTION
* ---------------------------------
*
* We implement a nested constant elasticity of substitution (CES) production function with capital, labor, and the (commercial) end-use services 
* represented in MESSAGE as inputs. :math:`\text{Y}_{n,t}` should correspond to gross domestic product (GDP).
*
* .. math:: \text{Y}_{n,t} = \left( a_{n} \cdot \text{K}_{n, t}^{ ( \rho_{n} \cdot \alpha_{n} ) } \cdot \text{L}_{n, t}^{ ( \rho_{n} \cdot ( 1 - \alpha_{n} ) ) } + \sum_{s} ( b_{n, s} \cdot \text{YE}_{n, s, t}^{\rho_{n}} ) \right)^{ \frac{1}{\rho_{n}} } \qquad \forall n, t > 1 \qquad (16)
*
***

PRODUCTION(node_active, year) $ (NOT macro_base_period(year))..
Y(node_active, year) =E=
    (LAKL(node_active) * K(node_active, year)**(rho(node_active) * kpvs(node_active)) * labor(node_active, year)**(rho(node_active) * (1 - kpvs(node_active)))
    + SUM(sector, PRFCONST(node_active, sector) * YE(node_active, sector, year)**rho(node_active))
    )**(1/rho(node_active))
;

$ontext
TECHNOLOGY(node_active, year) $ (NOT macro_base_period(year))..
AE(node_active, year) =E=
(gdp_calibrate(node_active, year)/ 1000)
/ ( (LAKL(node_active) * K(node_active, year)**(rho(node_active) * kpvs(node_active)) * labor(node_active, year)**(rho(node_active) * (1 - kpvs(node_active)))
+ PRFCONST(node_active, 'i_spec') * YE(node_active, 'i_spec', year)**rho(node_active)
+ PRFCONST(node_active, 'i_therm') * YE(node_active, 'i_therm', year)**rho(node_active)
+ PRFCONST(node_active, 'rc_spec') * YE(node_active, 'rc_spec', year)**rho(node_active)
+ PRFCONST(node_active, 'rc_therm') * YE(node_active, 'rc_therm', year)**rho(node_active)
+ PRFCONST(node_active, 'transport') * YE(node_active, 'transport', year)**rho(node_active)
)**(1/rho(node_active)) )
;
$offtext

MARGINAL_PRODUCT_CAPITAL(node_active, year) $ (NOT macro_base_period(year))..
INTEREST(node_active, year) =E= 
    Y(node_active, year)**(1-rho(node_active)) * LAKL(node_active) * kpvs(node_active) 
    * K(node_active, year)**(rho(node_active)*kpvs(node_active) - 1) 
    * labor(node_active, year)**(rho(node_active)*(1-kpvs(node_active)))
;

MARGINAL_PRODUCT_LABOR(node_active, year) $ (NOT macro_base_period(year))..
WAGE(node_active, year) =E=
    Y(node_active, year)**(1-rho(node_active)) * LAKL(node_active) * (1-kpvs(node_active)) 
    * K(node_active, year)**(rho(node_active)*kpvs(node_active)) 
    * labor(node_active, year)**(rho(node_active)*(1-kpvs(node_active)) - 1)
;

$ontext
MARGINAL_PRODUCT_CAPITAL(node_active, year) $ (NOT macro_base_period(year))..
INTEREST(node_active, year) + depr(node_active) =E=
Y(node_active, year) * LAKL(node_active) * kpvs(node_active)
* K(node_active, year)**(rho(node_active)*kpvs(node_active) - 1) 
* labor(node_active, year)**(rho(node_active)*(1 - kpvs(node_active)))
/ (LAKL(node_active) * K(node_active, year)**(rho(node_active) * kpvs(node_active)) 
    * labor(node_active, year)**(rho(node_active) * (1 - kpvs(node_active)))
    + PRFCONST(node_active, 'i_spec') * YE(node_active, 'i_spec', year)**rho(node_active)
    + PRFCONST(node_active, 'i_therm') * YE(node_active, 'i_therm', year)**rho(node_active)
    + PRFCONST(node_active, 'rc_spec') * YE(node_active, 'rc_spec', year)**rho(node_active)
    + PRFCONST(node_active, 'rc_therm') * YE(node_active, 'rc_therm', year)**rho(node_active)
    + PRFCONST(node_active, 'transport') * YE(node_active, 'transport', year)**rho(node_active))
;

MARGINAL_PRODUCT_LABOR(node_active, year) $ (NOT macro_base_period(year))..
WAGE(node_active, year) =E=
Y(node_active, year) * LAKL(node_active) * (1-kpvs(node_active))
* K(node_active, year)**(rho(node_active)*kpvs(node_active)) 
* labor(node_active, year)**(rho(node_active)*(1-kpvs(node_active))-1)
/ (LAKL(node_active) * K(node_active, year)**(rho(node_active) * kpvs(node_active)) 
    * labor(node_active, year)**(rho(node_active) * (1 - kpvs(node_active)))
    + PRFCONST(node_active, 'i_spec') * YE(node_active, 'i_spec', year)**rho(node_active)
    + PRFCONST(node_active, 'i_therm') * YE(node_active, 'i_therm', year)**rho(node_active)
    + PRFCONST(node_active, 'rc_spec') * YE(node_active, 'rc_spec', year)**rho(node_active)
    + PRFCONST(node_active, 'rc_therm') * YE(node_active, 'rc_therm', year)**rho(node_active)
    + PRFCONST(node_active, 'transport') * YE(node_active, 'transport', year)**rho(node_active))
;
$offtext

***
* Equation CAPITAL_GROWTH
* ---------------------------------
* Definition of the capital growth rate KGROW as the net growth rate of capital stock (backwards-looking).
* This is used in the terminal condition to ensure sufficient investment for capital replacement and growth.
* The growth rate is defined as: KGROW(t) = (K(t) - K(t-1)) / K(t-1)
*
* .. math:: \text{KGROW}_{n,t} = \frac{K_{n,t} - K_{n,t-1}}{K_{n,t-1}}
*
***

CAPITAL_GROWTH(node_active, year) $ (NOT macro_base_period(year))..
KGROW(node_active, year) =E=
(K(node_active, year) - SUM(year2$( seq_period(year2, year) ), K(node_active, year2))) / SUM(year2$( seq_period(year2, year) ), K(node_active, year2))
;

***
* Equations ENERGY ACCOUNTING
* ---------------------------------
* Energy accounting equations. See model documentation.
*
***

ENERGY_ACCOUNTING(node_active, sector, year) $ (NOT macro_base_period(year))..
TE(node_active, sector, year) =E=
YE(node_active, sector, year) + E(node_active, sector, year)
;

ENERGY_ACCOUNTING2(node_active, sector, year) $ (NOT macro_base_period(year))..
E(node_active, sector, year) =E=
TE(node_active, sector, year) * h(node_active, sector)
;

$ontext
TE_EQUATION(node_active, sector, year) $ (NOT macro_base_period(year))..
E(node_active, sector, year) =E=

( (labor(node_active, year) * EMIN(node_active) 
 + (beta_rc_spec(node_active) / alpha(node_active)) * C(node_active, year) / (eneprice(node_active, 'rc_spec', year)/1000)
) / h(node_active, 'rc_spec') ) $ (sameas(sector,'rc_spec') AND h(node_active, 'rc_spec') <> 0)

+ ( (labor(node_active, year) * EMIN(node_active)
   + (beta_rc_therm(node_active) / alpha(node_active)) * C(node_active, year) / (eneprice(node_active, 'rc_therm', year)/1000)
) / h(node_active, 'rc_therm') ) $ (sameas(sector,'rc_therm') AND h(node_active, 'rc_therm') <> 0)

+ ( (labor(node_active, year) * EMIN(node_active)
   + (beta_transport(node_active) / alpha(node_active)) * C(node_active, year) / (eneprice(node_active, 'transport', year)/1000)
) / h(node_active, 'transport') ) $ (sameas(sector,'transport') AND h(node_active, 'transport') <> 0)

+ 0 $ (sameas(sector,'i_spec') OR sameas(sector,'i_therm'))
;
$offtext

ENERGY_SUPPLY(node_active, sector, year) $ (NOT macro_base_period(year))..
PHYSENE(node_active, sector, year) =G=
TE(node_active, sector, year) * aeei_factor(node_active, sector, year)
;

COST_ENERGY(node_active, year) $ (NOT macro_base_period(year))..
EC(node_active, year) =E=
(total_cost(node_active, year)/1000
+ SUM(sector, eneprice(node_active, sector, year) * 1E-3 * (PHYSENE(node_active, sector, year) - enestart(node_active, sector, year)))
+ SUM(sector, eneprice(node_active, sector, year) * 1E-3 / enestart(node_active, sector, year) * (PHYSENE(node_active, sector, year) - enestart(node_active, sector, year)) * (PHYSENE(node_active, sector, year) - enestart(node_active, sector, year))))
;

***
* Equation TERMINAL_CONDITION
* ---------------------------------
* Given the finite time horizon of MACRO, a terminal constraint needs to be applied to ensure that investments are chosen at an appropriate level, i.e. to replace depriciated capital and
* provide net growth of capital stock beyond MACRO's time horizon :cite:`Manne-Richels-1992`. The goal is to avoid to the extend possible model artifacts resulting from this finite time horizon
* cutoff.
*
* .. math:: \text{K}_{n, y} \cdot  \left( \text{grow}_{n, y} + \text{depr}_n \right) \leq \text{I}_{n, y} \qquad \forall{ n, y = \text{last year}}
***

TERMINAL_CONDITION(node_active, last_period)..
I(node_active, last_period) =G= K(node_active, last_period) * (KGROW(node_active, last_period) + depr(node_active))
;

* ------------------------------------------------------------------------------
* model definition
* ------------------------------------------------------------------------------

MODEL MESSAGE_MACRO /
    UTILITY_FUNCTION
    CAPITAL_CONSTRAINT
    CAPITAL
    INVESTMENT
    PRODUCTION
*TE_EQUATION
    ENERGY_ACCOUNTING
    ENERGY_ACCOUNTING2
    ENERGY_SUPPLY
    COST_ENERGY
    TERMINAL_CONDITION
    MARGINAL_PRODUCT_CAPITAL
    MARGINAL_PRODUCT_LABOR
    CAPITAL_GROWTH
*TECHNOLOGY
/ ;

MESSAGE_MACRO.optfile = 1;
