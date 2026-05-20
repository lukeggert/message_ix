
* ------------------------------------------------------------------------------
* Start values of variables inside feasible domain (positive variables)
* ------------------------------------------------------------------------------

SVKN(node_macro, year) = ((potential_gdp(node_macro, year) - SUM(year2$( seq_period(year2,year) ), potential_gdp(node_macro, year2)) * depr(node_macro)**duration_period(year)) * kgdp(node_macro)) $ (NOT macro_base_period(year));

* NB This value is an *estimate* provided to the solver as a initial value (NEWENE.L in the next statement) in the
* problem domain. Depending on the interpretation of NEWENE (for instance, as an instantaneous value at the *start* or
* *end* of a period, or some point between; or a mean/average over the whole period or the representative year) and of
* base_demand, the value may be low or high. This affects solver performance but should not affect the optimal solution.
* See also:
* - The documentation of duration_period_sum.
* - Comments at https://github.com/iiasa/message_ix/pull/926 and the related issue #925.

SVNEWE(node_macro, sector, year) = (
  demand_base(node_macro, sector) * growth_factor(node_macro, year)
  - demand_base(node_macro, sector) * depr(node_macro) ** (
    SUM(year2$macro_base_period(year2), duration_period_sum(year2, year) + duration_period(year))
  )
)$(NOT macro_base_period(year));

*TE.L(node_macro, sector, macro_horizon) = (
*  SVNEWE(node_macro, sector, macro_horizon)$(SVNEWE(node_macro, sector, macro_horizon) > 0) + epsilon
*);
*YE.L(node_macro, sector, macro_horizon) = (1-h(node_macro, sector)) * (
*  SVNEWE(node_macro, sector, macro_horizon)$(SVNEWE(node_macro, sector, macro_horizon) > 0) + epsilon
*);
*E.L(node_macro, sector, macro_horizon) = h(node_macro, sector) * (
*  SVNEWE(node_macro, sector, macro_horizon)$(SVNEWE(node_macro, sector, macro_horizon) > 0) + epsilon
*);

PHYSENE.L(node_macro, sector, year)  = enestart(node_macro, sector, year) ;
K.L(node_macro, macro_horizon)  = SVKN(node_macro, macro_horizon) $ (SVKN(node_macro, macro_horizon) > 0) + epsilon ;

KGROW.L(node_macro, macro_base_period) = grow(node_macro, macro_base_period) ; 

* ------------------------------------------------------------------------------
* Lower bounds on variables help to avoid singularities
* ------------------------------------------------------------------------------

K.LO(node_macro, macro_horizon)  = LOTOL(node_macro) * k0(node_macro) ;
Y.LO(node_macro, macro_horizon)  = LOTOL(node_macro) * y0(node_macro) ;

C.LO(node_macro, macro_horizon)  = LOTOL(node_macro) * c0(node_macro) ;
I.LO(node_macro, macro_horizon)  = LOTOL(node_macro) * i0(node_macro) ;

TE.LO(node_macro, sector, macro_horizon) = LOTOL(node_macro) * enestart(node_macro, sector, macro_horizon) / aeei_factor(node_macro, sector, macro_horizon) ;
YE.LO(node_macro, sector, macro_horizon) = (1-h(node_macro, sector)) * LOTOL(node_macro) * enestart(node_macro, sector, macro_horizon) / aeei_factor(node_macro, sector, macro_horizon) ;
E.LO(node_macro, sector, macro_horizon) = h(node_macro, sector) * LOTOL(node_macro) * enestart(node_macro, sector, macro_horizon) / aeei_factor(node_macro, sector, macro_horizon) ;

* Sicherheitsprüfung: Korrektur zu hoher Lower Bounds (verhindert inconsistent bounds error)
*TE.LO(node_macro, sector, macro_horizon)$(TE.LO(node_macro, sector, macro_horizon) > 1e6) = 1e-6;
*YE.LO(node_macro, sector, macro_horizon)$(YE.LO(node_macro, sector, macro_horizon) > 1e6) = 1e-6;
*E.LO(node_macro, sector, macro_horizon)$(E.LO(node_macro, sector, macro_horizon) > 1e6) = 1e-6;

AE.LO(node_macro, macro_horizon) = 0.1 ; 
WAGE.LO(node_macro, macro_horizon) = 0.01 ; 
*INTEREST.LO(node_macro, macro_horizon) = 0.02 ; 

* ------------------------------------------------------------------------------
* Upper bounds on variables help to avoid singularities
* ------------------------------------------------------------------------------



* ------------------------------------------------------------------------------
* Base year values of variables are fixed to historical values
* ------------------------------------------------------------------------------

* division by aeei_factor is necesary in case MACRO starts after initialize_period (in case of slicing)
TE.FX(node_macro, sector, macro_base_period) = demand_base(node_macro, sector) / aeei_factor(node_macro, sector, macro_base_period) ;
YE.FX(node_macro, sector, macro_base_period) = (1-h(node_macro, sector)) * demand_base(node_macro, sector) / aeei_factor(node_macro, sector, macro_base_period) ;
E.FX(node_macro, sector, macro_base_period) = h(node_macro, sector) * demand_base(node_macro, sector) / aeei_factor(node_macro, sector, macro_base_period) ;

Y.FX(node_macro, macro_base_period) = y0(node_macro) ;
K.FX(node_macro, macro_base_period) = k0(node_macro) ;
C.FX(node_macro, macro_base_period) = c0(node_macro) ;
I.FX(node_macro, macro_base_period) = i0(node_macro) ;
EC.FX(node_macro, macro_base_period) = y0(node_macro) - i0(node_macro) - c0(node_macro) ; 

AE.FX(node_macro, macro_base_period) = ( gdp_base(node_macro) / (lakl(node_macro) * k0(node_macro)**(rho(node_macro) * kpvs(node_macro)) * 
  labor(node_macro, macro_base_period)**(rho(node_macro) * (1 - kpvs(node_macro))) +
  PRFCONST(node_macro, 'i_spec') * YE.L(node_macro, 'i_spec', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'i_therm') * YE.L(node_macro, 'i_therm', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'rc_spec') * YE.L(node_macro, 'rc_spec', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'rc_therm') * YE.L(node_macro, 'rc_therm', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'transport') * YE.L(node_macro, 'transport', macro_base_period)**rho(node_macro))**(1/rho(node_macro)) ) 
;

INTEREST.FX(node_macro, macro_base_period) = y0(node_macro) * lakl(node_macro) *
  kpvs(node_macro) * k0(node_macro)**(rho(node_macro)*kpvs(node_macro) -1) *
  labor(node_macro, macro_base_period)**(rho(node_macro)*(1-kpvs(node_macro))) * 
  (lakl(node_macro) * k0(node_macro)**(rho(node_macro)*kpvs(node_macro)) *
  labor(node_macro, macro_base_period)**(rho(node_macro)*(1-kpvs(node_macro))) +
  PRFCONST(node_macro, 'i_spec') * YE.L(node_macro, 'i_spec', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'i_therm') * YE.L(node_macro, 'i_therm', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'rc_spec') * YE.L(node_macro, 'rc_spec', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'rc_therm') * YE.L(node_macro, 'rc_therm', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'transport') * YE.L(node_macro, 'transport', macro_base_period)**rho(node_macro))**(-1) - depr(node_macro)
;

WAGE.FX(node_macro, macro_base_period) =  y0(node_macro) * lakl(node_macro) *
  (1-kpvs(node_macro)) * k0(node_macro)**(rho(node_macro)*kpvs(node_macro)) *
  labor(node_macro, macro_base_period)**(rho(node_macro)*(1-kpvs(node_macro))-1) * 
  (lakl(node_macro) * k0(node_macro)**(rho(node_macro)*kpvs(node_macro)) *
  labor(node_macro, macro_base_period)**(rho(node_macro)*(1-kpvs(node_macro))) +
  PRFCONST(node_macro, 'i_spec') * YE.L(node_macro, 'i_spec', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'i_therm') * YE.L(node_macro, 'i_therm', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'rc_spec') * YE.L(node_macro, 'rc_spec', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'rc_therm') * YE.L(node_macro, 'rc_therm', macro_base_period)**rho(node_macro) +
  PRFCONST(node_macro, 'transport') * YE.L(node_macro, 'transport', macro_base_period)**rho(node_macro))**(-1) 
;

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
