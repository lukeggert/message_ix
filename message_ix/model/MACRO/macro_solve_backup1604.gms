
* ------------------------------------------------------------------------------
* Start values of variables inside feasible domain (positive variables)
* ------------------------------------------------------------------------------

*SVKN(node_macro, year) = ((gdp_calibrate(node_macro, year)/1000) - SUM(year2$( seq_period(year2,year) ), (gdp_calibrate(node_macro, year2)/1000)) * depr(node_macro)**duration_period(year) * kgdp(node_macro)) $ (NOT macro_base_period(year));
*SVKN(node_macro, year) = k0(node_macro) * (1+0.05)**duration_period(year) ; 
SVKN(node_macro, year) = kgdp(node_macro) * gdp_calibrate(node_macro, year) / 1000; 

* NB This value is an *estimate* provided to the solver as a initial value (NEWENE.L in the next statement) in the
* problem domain. Depending on the interpretation of NEWENE (for instance, as an instantaneous value at the *start* or
* *end* of a period, or some point between; or a mean/average over the whole period or the representative year) and of
* base_demand, the value may be low or high. This affects solver performance but should not affect the optimal solution.
* See also:
* - The documentation of duration_period_sum.
* - Comments at https://github.com/iiasa/message_ix/pull/926 and the related issue #925.

SVTE(node_macro, sector, year) = (
  demand_base(node_macro, sector) * growth_factor(node_macro, year)
  - demand_base(node_macro, sector) * depr(node_macro) ** (
    SUM(year2$macro_base_period(year2), duration_period_sum(year2, year) + duration_period(year))
  )
)$(NOT macro_base_period(year));

SVNEWE(node_macro, sector, year) = (
  demand_base(node_macro, sector) * growth_factor(node_macro, year)
  - demand_base(node_macro, sector) * (1 - depr(node_macro)) ** (
    SUM(year2$macro_base_period(year2), duration_period_sum(year2, year) + duration_period(year))
  )
)$(NOT macro_base_period(year));

NEWENE.L(node_macro, sector, macro_horizon) = (
  SVNEWE(node_macro, sector, macro_horizon)$(SVNEWE(node_macro, sector, macro_horizon) > 0) + epsilon
);

TE.L(node_macro, sector, macro_horizon) = (
  SVTE(node_macro, sector, macro_horizon)$(SVTE(node_macro, sector, macro_horizon) > 0) + epsilon
);
YE.L(node_macro, sector, macro_horizon) = (1-h(node_macro, sector)) * (
  SVTE(node_macro, sector, macro_horizon)$(SVTE(node_macro, sector, macro_horizon) > 0) + epsilon
);
E.L(node_macro, sector, macro_horizon) = h(node_macro, sector) * (
  SVTE(node_macro, sector, macro_horizon)$(SVTE(node_macro, sector, macro_horizon) > 0) + epsilon
);

*KGROW.L(node_macro, macro_horizon) = grow(node_macro, macro_horizon) ; 

*C.L(node_macro, macro_horizon) = c0(node_macro) ; 

I.L(node_macro, macro_horizon) = SVKN(node_macro, macro_horizon) * (grow(node_macro, macro_horizon) + depr(node_macro)) ; 
K.L(node_macro, macro_horizon)  = SVKN(node_macro, macro_horizon) $ (SVKN(node_macro, macro_horizon) > 0) + epsilon ; 
PHYSENE.L(node_macro, sector, year)  = enestart(node_macro, sector, year) ;
TE.L(node_macro, sector, year) = enestart(node_macro, sector, year) / aeei_factor(node_macro, sector, year);
YE.L(node_macro, sector, year) = (1-h(node_macro, sector)) * enestart(node_macro, sector, year) / aeei_factor(node_macro, sector, year);
E.L(node_macro, sector, year) = h(node_macro, sector) * enestart(node_macro, sector, year) / aeei_factor(node_macro, sector, year);

C.L(node_macro, macro_horizon) = gdp_calibrate(node_macro, macro_horizon)/ 1000 - (SVKN(node_macro, macro_horizon) * (grow(node_macro, macro_horizon) + depr(node_macro))) - ecst0(node_macro)/1000 ;
Y.L(node_macro, macro_horizon) = gdp_calibrate(node_macro, macro_horizon) / 1000 ; 

$ontext
E.L(node_macro, sector, macro_horizon) $ sameas(sector, 'rc_spec') =
  labor(node_macro, macro_horizon) * EMIN(node_macro)
  + (beta_rc_spec(node_macro) / alpha(node_macro)) * C.L(node_macro, macro_horizon) / (eneprice(node_macro, 'rc_spec', macro_horizon)/1000) ;
E.L(node_macro, sector, macro_horizon) $ sameas(sector, 'rc_therm') =
  labor(node_macro, macro_horizon) * EMIN(node_macro)
  + (beta_rc_therm(node_macro) / alpha(node_macro)) * C.L(node_macro, macro_horizon) / (eneprice(node_macro, 'rc_therm', macro_horizon)/1000) ;
E.L(node_macro, sector, macro_horizon) $ sameas(sector, 'transport') =
  labor(node_macro, macro_horizon) * EMIN(node_macro)
  + (beta_transport(node_macro) / alpha(node_macro)) * C.L(node_macro, macro_horizon) / (eneprice(node_macro, 'transport', macro_horizon)/1000) ;
E.L(node_macro, sector, macro_horizon) $ (sameas(sector, 'i_spec') OR sameas(sector, 'i_therm')) = 0 ;
YE.L(node_macro, sector, macro_horizon) = MAX(epsilon, TE.L(node_macro, sector, macro_horizon) - E.L(node_macro, sector, macro_horizon)) ;
$offtext

INTEREST.L(node_macro, macro_horizon) = 0.1 ; 
WAGE.L(node_macro, macro_horizon) = 1 ; 

* ------------------------------------------------------------------------------
* Lower bounds on variables help to avoid singularities
* ------------------------------------------------------------------------------

K.LO(node_macro, macro_horizon)  = k0(node_macro) ;
Y.LO(node_macro, macro_horizon)  = LOTOL(node_macro) * y0(node_macro) ;

C.LO(node_macro, macro_horizon)  = LOTOL(node_macro) * c0(node_macro) ;

I.LO(node_macro, macro_horizon)  = LOTOL(node_macro) * i0(node_macro) ;

TE.LO(node_macro, sector, macro_horizon) = LOTOL(node_macro) * enestart(node_macro, sector, macro_horizon) ;
YE.LO(node_macro, sector, macro_horizon) = MAX(
  (1-h(node_macro, sector)) * enestart(node_macro, sector, macro_horizon),
  SUM(year2$seq_period(year2, macro_horizon), (1-h(node_macro, sector)) * enestart(node_macro, sector, year2)) + epsilon
) ;
E.LO(node_macro, sector, macro_horizon) = h(node_macro, sector) * LOTOL(node_macro) * enestart(node_macro, sector, macro_horizon) ;

NEWENE.LO(node_macro, sector, macro_horizon) = 0.001 ; 

*TE.LO(node_macro, sector, macro_horizon) = demand_base(node_macro, sector) / aeei_factor(node_macro, sector, macro_horizon) ;
*YE.LO(node_macro, sector, macro_horizon) = (1-h(node_macro, sector)) * demand_base(node_macro, sector) / aeei_factor(node_macro, sector, macro_horizon) ;
*E.LO(node_macro, sector, macro_horizon) = h(node_macro, sector) * demand_base(node_macro, sector) / aeei_factor(node_macro, sector, macro_horizon) ;

*AE.LO(node_macro, macro_horizon) = 0.1 ; 
WAGE.LO(node_macro, macro_horizon) = 0.01 ; 
INTEREST.LO(node_macro, macro_horizon) = 0.03 ; 

*EC.LO(node_macro, macro_horizon) = (y0(node_macro) - i0(node_macro) - c0(node_macro)) ;

* ------------------------------------------------------------------------------
* Upper bounds
* ------------------------------------------------------------------------------

*PHYSENE.UP(node_macro, sector, year) = enestart(node_macro, sector, year) ;

* Scale NEWENE cap with model time-step length (5y -> 1.0, 10y -> 2.0)
newene_share_cap(node_macro, sector, year) = duration_period(year) / 7 ; 

* / 7 oben für R12_AFR; alpha 0.82, betas 0.06 (oder ohne 7 oben; alpha 0.88, beta_transport = 0.06, betas = 0.03)
* ohne 7 oben; alpha 0.88, beta_transport = 0.06, betas = 0.03

*NEWENE.UP(node_macro, sector, year) $ (NOT macro_base_period(year)) =
*    min(
*      newene_share_cap(node_macro, sector, year)
*      * SUM(year2$(seq_period(year2,year)), TE.L(node_macro, sector, year2)),
*      3 * newene_share_cap(node_macro, sector, year)
*      * SUM(year2$(seq_period(year2,year)), PHYSENE.L(node_macro, sector, year2))
*    ) ;

*NEWENE.UP(node_macro, sector, year) $ (NOT macro_base_period(year)) =
*      newene_share_cap(node_macro, sector, year)
*      * SUM(year2$(seq_period(year2,year)), TE.L(node_macro, sector, year2))
*      ;

*TE.UP(node_macro, sector, year) $ (NOT macro_base_period(year)) =
*  1.5 * enestart(node_macro, sector, year) ;

* ------------------------------------------------------------------------------
* Base year values of variables are fixed to historical values
* ------------------------------------------------------------------------------

* qunatile values

* When different betas and alphas: ohne lower für cap und con optimal solution, aber single values for a specific quantile. With lower no solution
KAP.LO(node_macro, macro_horizon, quantile) = LOTOL(node_macro) * k0(node_macro) * quantile_share(quantile) ;
CON.LO(node_macro, macro_horizon, quantile) = LOTOL(node_macro) * c0(node_macro) * quantile_share(quantile) ;
*LAB.LO(node_macro, macro_horizon, quantile) = 0 ;


*KAP.FX(node_macro, macro_base_period, quantile) = k0(node_macro) * quantile_share(quantile) ;
*CON.FX(node_macro, macro_base_period, quantile) = c0(node_macro) * quantile_share(quantile) ;
*LAB.FX(node_macro, macro_base_period, quantile) = labor(node_macro, macro_base_period) * quantile_share(quantile) ;

KAP.L(node_macro, macro_horizon, quantile) = K.L(node_macro, macro_horizon) * quantile_share(quantile) ;
CON.L(node_macro, macro_horizon, quantile) = C.L(node_macro, macro_horizon) * quantile_share(quantile) ;
LAB.L(node_macro, macro_base_period, quantile) = labor(node_macro, macro_base_period) * quantile_share(quantile) ;

* satrtwerte mit quantilfsverteilung funktioniert. FX auf 2025 setzen führt zu infeasibilty.

* DLE:

*E_Q.LO(node_macro, sector, macro_horizon, quantile) $
*  (sameas(sector, 'rc_spec') OR sameas(sector, 'rc_therm') OR sameas(sector, 'transport')) =
*  e_q_min(node_macro, sector, macro_horizon) ;

*E_Q.L(node_macro, sector, year, quantile) $
*    (sameas(sector, 'rc_spec') OR sameas(sector, 'rc_therm') OR sameas(sector, 'transport')) = max(
*  e_q_min(node_macro, sector, year),
*  E.L(node_macro, sector, year) * quantile_share(quantile)
*) ;


* division by aeei_factor is necesary in case MACRO starts after initialize_period (in case of slicing)
TE.FX(node_macro, sector, macro_base_period) = demand_base(node_macro, sector) / aeei_factor(node_macro, sector, macro_base_period) ;
YE.FX(node_macro, sector, macro_base_period) = (1-h(node_macro, sector)) * demand_base(node_macro, sector) / aeei_factor(node_macro, sector, macro_base_period) ;
E.FX(node_macro, sector, macro_base_period) = h(node_macro, sector) * demand_base(node_macro, sector) / aeei_factor(node_macro, sector, macro_base_period) ;

Y.FX(node_macro, macro_base_period) = y0(node_macro) ;
K.FX(node_macro, macro_base_period) = k0(node_macro) ;
*EC.FX(node_macro, macro_base_period) = y0(node_macro) - i0(node_macro) - c0(node_macro) ;
EC.FX(node_macro, macro_base_period) = ecst0(node_macro) / 1000 ; 
I.FX(node_macro, macro_base_period) = i0(node_macro) ;
C.FX(node_macro,macro_base_period) = c0(node_macro) ; 

$ontext
AE.FX(node_macro, macro_base_period) = ( gdp_base(node_macro) / (lakl(node_macro) * k0(node_macro)**(rho(node_macro) * kpvs(node_macro)) * 
  labor(node_macro, macro_base_period)**(rho(node_macro) * (1 - kpvs(node_macro))) +
  PRFCONST(node_macro, 'i_spec') * YE.L(node_macro, 'i_spec', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'i_therm') * YE.L(node_macro, 'i_therm', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'rc_spec') * YE.L(node_macro, 'rc_spec', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'rc_therm') * YE.L(node_macro, 'rc_therm', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'transport') * YE.L(node_macro, 'transport', macro_base_period)**rho(node_macro))**(1/rho(node_macro)) ) 
;
$offtext

INTEREST.FX(node_macro, macro_base_period) = y0(node_macro)**(1-rho(node_macro)) * lakl(node_macro) * kpvs(node_macro) * 
  k0(node_macro)**(rho(node_macro)*kpvs(node_macro) -1) *
  labor(node_macro, macro_base_period)**(rho(node_macro)*(1-kpvs(node_macro)))
;

WAGE.FX(node_macro, macro_base_period) =  y0(node_macro)**(1-rho(node_macro)) * lakl(node_macro) * (1-kpvs(node_macro)) * 
  k0(node_macro)**(rho(node_macro)*kpvs(node_macro)) *
  labor(node_macro, macro_base_period)**(rho(node_macro)*(1-kpvs(node_macro))-1)
; 

$ontext
I.FX(node_macro, macro_base_period) = ((k0(node_macro) * ((1 + INTEREST.L(node_macro, macro_base_period))**duration_period(macro_base_period) - 1) / duration_period(macro_base_period))
  + labor(node_macro, macro_base_period) * WAGE.L(node_macro, macro_base_period)  
  - eneprice(node_macro, 'rc_spec', macro_base_period)/1000 * labor(node_macro, macro_base_period) * EMIN(node_macro)
  - eneprice(node_macro, 'rc_therm', macro_base_period)/1000 * labor(node_macro, macro_base_period) * EMIN(node_macro) 
  - eneprice(node_macro, 'transport', macro_base_period)/1000 * labor(node_macro, macro_base_period) * EMIN(node_macro) 
  - (((alpha(node_macro) + beta_rc_spec(node_macro) + beta_rc_therm(node_macro) + beta_transport(node_macro))/alpha(node_macro))) * (y0(node_macro) - EC.L(node_macro,macro_base_period))) 
  * (alpha(node_macro)/(alpha(node_macro) - (alpha(node_macro) + beta_rc_spec(node_macro) + beta_rc_therm(node_macro) + beta_transport(node_macro))))
;
C.FX(node_macro, macro_base_period) = y0(node_macro) - I.L(node_macro, macro_base_period) - EC.L(node_macro,macro_base_period) ; 
$offtext

$IFTHEN %MACRO_CONCURRENT% == "0"

DISPLAY "Solve MACRO for each node in sequence";

node_active(node) = NO ;

* Switch between solving for all or a specific region:

*LOOP(node$node_macro(node),
*  node_active(node_macro) = NO ;
*  node_active(node) = YES ;
*  DISPLAY node_active ;

LOOP(node$ (node_macro(node) AND sameas(node, "R12_AFR")),
  node_active(node_macro) = NO ;
  node_active(node) = YES ;
*  DISPLAY node_active ;

  SOLVE MESSAGE_MACRO MAXIMIZING UTILITY USING NLP ;

* Write model status summary for the current node
*  status(node,'modelstat') = MESSAGE_MACRO.modelstat ;
*  status(node,'solvestat') = MESSAGE_MACRO.solvestat ;
*  status(node,'resUsd')    = MESSAGE_MACRO.resUsd ;
*  status(node,'objEst')    = MESSAGE_MACRO.objEst ;
*  status(node,'objVal')    = MESSAGE_MACRO.objVal ;
);

$ELSE

DISPLAY "Solve MACRO for all nodes concurrently";

node_active(node_macro) = YES;

SOLVE MESSAGE_MACRO MAXIMIZING UTILITY USING NLP;

* Write model status summary for all nodes
* status('all','modelstat') = MESSAGE_MACRO.modelstat;
* status('all','solvestat') = MESSAGE_MACRO.solvestat;
* status('all','resUsd')    = MESSAGE_MACRO.resUsd;
* status('all','objEst')    = MESSAGE_MACRO.objEst;
* status('all','objVal')    = MESSAGE_MACRO.objVal;

$ENDIF
