$TITLE MACRO core MCP formulation with 10 Haushaltsdezi
$EOLCOM #

* Sets
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

* Parameter, Variable, Scalar, etc. Definitionen wie im Originalmodell
* ...

* --- MCP Variablen und Gleichungen ---

POSITIVE VARIABLES
    KAP(node_macro, year_all, quantile)     'Kapital pro Haushalt'
    CON(node_macro, year_all, quantile)     'Konsum pro Haushalt'
    I_HH(node_macro, year_all, quantile)    'Investition pro Haushalt'
    K(node_macro, year_all)                 'Aggregiertes Kapital'
    C(node_macro, year_all)                 'Aggregierter Konsum'
    I(node_macro, year_all)                 'Aggregierte Investition'
    Y(node_macro, year_all)                 'Produktion'
    WAGE(node_macro, year_all)              'Lohn'
    INTEREST(node_macro, year_all)          'Zins'
    YE(node_macro, sector, year_all)           'Produktiver Endnutzenergieinput'
    TE(node_macro, sector, year_all)           'Gesamter Endnutzenergieinput'
    E(node_macro, sector, year_all)            'Direkter Energiebedarf Haushalt'
    PHYSENE(node_macro, sector, year_all)      'Physische Endnutzenergie'
    EC(node_macro, year_all)                   'Systemkosten'
    mu_ENERGY_ACCOUNTING(node_macro, sector, year_all)
    mu_ENERGY_ACCOUNTING2(node_macro, sector, year_all)
    mu_ENERGY_SUPPLY(node_macro, sector, year_all)
    mu_COST_ENERGY(node_macro, year_all)      'Kosten für Energie'
    KGROW(node_macro, year_all)                  'Kapitalwachstum'
    HH_UTILITY_REP                               'Haushaltsnutzenfunktion (repräsentativ)'
        MU_CAPITAL(node_macro, year_all, quantile)      'Dual zu Kapitaldynamik Haushalt (FOC_KAP)'
        NU_TERM(node_macro, year_all, quantile)         'Dual zu Terminalbedingung Haushalt'
;

FREE VARIABLES
    MU_KAP(node_macro, year_all, quantile)  'Dual zu Kapitaldynamik Haushalt'
    MU_INV(node_macro, year_all, quantile)  'Dual zu Investitionsgleichung Haushalt'
    LAMBDA_CAP(node_macro, year_all)        'Dual zu aggregierter Kapitalverwendung'
     MU_PROD(node_macro, year_all)              'Dual zu Produktionsfunktion'
     MU_MPK(node_macro, year_all)               'Dual zu FOC Kapital'
     MU_MPL(node_macro, year_all)               'Dual zu FOC Arbeit'
    * weitere Dualvariablen nach Bedarf
    mu_ENERGY_ACCOUNTING(node_macro, sector, year_all)
    mu_ENERGY_ACCOUNTING2(node_macro, sector, year_all)
    mu_ENERGY_SUPPLY(node_macro, sector, year_all)
    MU_AGG_KAP(node_macro, year_all)           'Dual zu AGG_KAP (Kapitalaggregation)'
    MU_AGG_CON(node_macro, year_all)           'Dual zu AGG_CON (Konsumaggregation)'
    MU_AGG_INV(node_macro, year_all)           'Dual zu AGG_INV (Investitionsaggregation)'
    MU_CAP_USE(node_macro, year_all)           'Dual zu CAP_USE (Kapitalverwendung)'
    MU_TERM_COND(node_macro)                   'Dual zu TERMINAL_CONDITION_MCP (Terminalbedingung Kapital)'
    MU_HH_UTILITY_REP                          'Dual zu HH_UTILITY_REP_DEF (repräsentative Nutzenfunktion)'
    MU_FOC_CON(node_macro, year_all, quantile) 'Dual zu FOC_CON (FOC Konsum Haushalt)'
    MU_FOC_KAP(node_macro, year_all, quantile) 'Dual zu FOC_KAP (FOC Kapital Haushalt)'
    MU_LABOR_MARKET(node_macro, year_all)      'Dual zu LABOR_MARKET (Arbeitsmarkt)'
    MU_TERMINAL_HH(node_macro, year_all, quantile) 'Dual zu TERMINAL_HH_MCP (Terminalbedingung Haushalt)'
    MU_KGROW(node_macro, year_all)             'Dual zu CAPITAL_GROWTH_MCP (Kapitalwachstum)'
    mu_COST_ENERGY(node_macro, year_all)
;

* --- MCP Gleichungen ---

EQUATIONS
    KAP_DYN(node_macro, year_all, quantile)     'Kapitaldynamik Haushalt'
    INV_ACC(node_macro, year_all, quantile)     'Investitionsgleichung Haushalt'
    AGG_KAP(node_macro, year_all)               'Aggregiertes Kapital'
    AGG_CON(node_macro, year_all)               'Aggregierter Konsum'
    AGG_INV(node_macro, year_all)               'Aggregierte Investition'
    CAP_USE(node_macro, year_all)               'Aggregierte Kapitalverwendung (Y = C + I + EC)'
    PROD_FUNC(node_macro, year_all)                'Produktionsfunktion'
    FOC_KAP_PROD(node_macro, year_all)             'FOC Kapital (MPK)'
    FOC_LAB_PROD(node_macro, year_all)             'FOC Arbeit (MPL)'
    * weitere Gleichungen nach Bedarf
    ENERGY_ACCOUNTING_MCP(node_macro, sector, year_all)   'Energie-Bilanzgleichung (MCP)'
    ENERGY_ACCOUNTING2_MCP(node_macro, sector, year_all)  'Haushalts-Energiebedarf (MCP)'
    ENERGY_SUPPLY_MCP(node_macro, sector, year_all)       'Energieangebot (MCP)'
    COST_ENERGY_MCP(node_macro, year_all)                 'Systemkosten (MCP)'
    TERMINAL_CONDITION_MCP(node_macro)                    'Terminalbedingung für Kapital im letzten Jahr'
    HH_UTILITY_REP_DEF                                    'Definition der Haushaltsnutzenfunktion (repräsentativ)'
        FOC_CON(node_macro, year_all, quantile)              'FOC für Haushaltskonsum pro quantile'
        FOC_KAP(node_macro, year_all, quantile)              'FOC für Haushaltskapital pro quantile'
        CAPITAL_MARKET(node_macro, year_all)                 'Kapitalmarktgleichgewicht'
        LABOR_MARKET(node_macro, year_all)                   'Arbeitsmarktgleichgewicht'
            TERMINAL_HH_MCP(node_macro, year_all, quantile)      'Terminalbedingung für Haushalte (letztes Jahr)'
;

* Beispielhafte Gleichungen (exakt, keine Vereinfachung):

* Produktionsfunktion (wie im Original)
PROD_FUNC(node_macro, year)$(NOT macro_base_period(year))..
    Y(node_macro, year) =E=
        (ACONST(node_macro) * K(node_macro, year)**(rho(node_macro) * kpvs(node_macro))
        * labor(node_macro, year)**(rho(node_macro) * (1 - kpvs(node_macro)))
        + SUM(sector, BCONST(node_macro, sector) * YE(node_macro, sector, year)**rho(node_macro))
        )**(1/rho(node_macro));

* FOC Kapital (MPK = Zins)
FOC_KAP_PROD(node_macro, year)$(NOT macro_base_period(year))..
    INTEREST(node_macro, year) =E=
        Y(node_macro, year)**(1 - rho(node_macro)) * LAKL(node_macro) * kpvs(node_macro)
        * K(node_macro, year)**(rho(node_macro)*kpvs(node_macro) - 1)
        * labor(node_macro, year)**(rho(node_macro)*(1 - kpvs(node_macro)));

* FOC Arbeit (MPL = Lohn)
FOC_LAB_PROD(node_macro, year)$(NOT macro_base_period(year))..
    WAGE(node_macro, year) =E=
        Y(node_macro, year)**(1 - rho(node_macro)) * LAKL(node_macro) * (1 - kpvs(node_macro))
        * K(node_macro, year)**(rho(node_macro)*kpvs(node_macro))
        * labor(node_macro, year)**(rho(node_macro)*(1 - kpvs(node_macro)) - 1);

KAP_DYN(node_macro, year, quantile)$(NOT macro_base_period(year))..
    KAP(node_macro, year, quantile) =E=
        SUM(year2$(seq_period(year2,year)),
            KAP(node_macro, year2, quantile) * (1 - depr(node_macro))**duration_period(year)
            + duration_period(year) * I_HH(node_macro, year, quantile)
        );

INV_ACC(node_macro, year, quantile)$(NOT macro_base_period(year) AND NOT last_period(year))..
    I_HH(node_macro, year, quantile) =E=
        SUM(year2$(seq_period(year2,year)),
            (KAP(node_macro, year2, quantile) * ((1 + INTEREST(node_macro, year))**duration_period(year) - 1) / duration_period(year))
            + LAB(node_macro, year, quantile) * WAGE(node_macro, year)
            - eneprice(node_macro, 'rc_spec', year)/1000 * quantile_share(quantile) * EMIN(node_macro)
            - eneprice(node_macro, 'rc_therm', year)/1000 * quantile_share(quantile) * EMIN(node_macro)
            - eneprice(node_macro, 'transport', year)/1000 * quantile_share(quantile) * EMIN(node_macro)
            - ((alpha_q(node_macro) + beta_rc_spec_q(node_macro) + beta_rc_therm_q(node_macro) + beta_transport_q(node_macro))/alpha_q(node_macro)) * CON(node_macro, year, quantile)
        );

AGG_KAP(node_macro, year)$(NOT macro_base_period(year))..
    K(node_macro, year) =E= SUM(quantile, KAP(node_macro, year, quantile));

AGG_CON(node_macro, year)$(NOT macro_base_period(year))..
    C(node_macro, year) =E= SUM(quantile, CON(node_macro, year, quantile));

AGG_INV(node_macro, year)$(NOT macro_base_period(year))..
    I(node_macro, year) =E= SUM(quantile, I_HH(node_macro, year, quantile));

* Aggregierte Kapitalverwendung (z.B. Y = C + I + EC)
CAP_USE(node_macro, year)$(NOT macro_base_period(year))..
    Y(node_macro, year) =E= C(node_macro, year) + I(node_macro, year) + EC(node_macro, year);

* --- Energiegleichungen (MCP-Form) ---

* Energie-Bilanzgleichung
ENERGY_ACCOUNTING_MCP(node_macro, sector, year_all)$(NOT macro_base_period(year_all))..
    TE(node_macro, sector, year_all) =E= YE(node_macro, sector, year_all) + E(node_macro, sector, year_all);

* Haushalts-Energiebedarf
ENERGY_ACCOUNTING2_MCP(node_macro, sector, year_all)$(NOT macro_base_period(year_all))..
    E(node_macro, sector, year_all) =E=
        ( EMIN(node_macro)
        + (beta_rc_spec(node_macro) / alpha(node_macro)) * C(node_macro, year_all) / (eneprice(node_macro, 'rc_spec', year_all)/1000)
        ) $ sameas(sector, 'rc_spec')
    + ( EMIN(node_macro)
        + (beta_rc_therm(node_macro) / alpha(node_macro)) * C(node_macro, year_all) / (eneprice(node_macro, 'rc_therm', year_all)/1000)
        ) $ sameas(sector, 'rc_therm')
    + ( EMIN(node_macro)
        + (beta_transport(node_macro) / alpha(node_macro)) * C(node_macro, year_all) / (eneprice(node_macro, 'transport', year_all)/1000)
        ) $ sameas(sector, 'transport')
    + 0 $ (sameas(sector, 'i_spec') OR sameas(sector, 'i_therm'));

* Energieangebot
ENERGY_SUPPLY_MCP(node_macro, sector, year_all)$(NOT macro_base_period(year_all))..
    PHYSENE(node_macro, sector, year_all) =G= TE(node_macro, sector, year_all) * aeei_factor(node_macro, sector, year_all);

* Systemkosten
COST_ENERGY_MCP(node_macro, year_all)$(NOT macro_base_period(year_all))..
    EC(node_macro, year_all) =E=
        (total_cost(node_macro, year_all)/1000
        + SUM(sector, eneprice(node_macro, sector, year_all) * 1E-3 * (PHYSENE(node_macro, sector, year_all) - enestart(node_macro, sector, year_all)))
        + SUM(sector, eneprice(node_macro, sector, year_all) * 1E-3 / enestart(node_macro, sector, year_all)
            * (PHYSENE(node_macro, sector, year_all) - enestart(node_macro, sector, year_all)) * (PHYSENE(node_macro, sector, year_all) - enestart(node_macro, sector, year_all)))
        );

* Terminalbedingung (letztes Jahr)
TERMINAL_CONDITION_MCP(node_macro)$(last_period(year_all))..
    I(node_macro, year_all) =E= SUM(year2$seq_period(year2, year_all), K(node_macro, year_all) * (KGROW(node_macro, year_all) + depr(node_macro)));

* Haushaltsnutzenfunktion (repräsentativ)
HH_UTILITY_REP_DEF..
    HH_UTILITY_REP =E=
    SUM(node_macro,
        hh_scale * (
            SUM((year_all, quantile)$(macro_horizon(year_all) AND NOT macro_base_period(year_all) AND NOT last_period(year_all)),
                udf(node_macro, year_all) * (
                    LOG(CON(node_macro, year_all, quantile))
                    - beta_rc_spec_q(node_macro, quantile) * LOG(eneprice(node_macro, 'rc_spec', year_all) / 1000)
                    + beta_rc_spec_q(node_macro, quantile) * LOG(beta_rc_spec_q(node_macro, quantile) / alpha_q(node_macro, quantile))
                    - beta_rc_therm_q(node_macro, quantile) * LOG(eneprice(node_macro, 'rc_therm', year_all) / 1000)
                    + beta_rc_therm_q(node_macro, quantile) * LOG(beta_rc_therm_q(node_macro, quantile) / alpha_q(node_macro, quantile))
                    - beta_transport_q(node_macro, quantile) * LOG(eneprice(node_macro, 'transport', year_all) / 1000)
                    + beta_transport_q(node_macro, quantile) * LOG(beta_transport_q(node_macro, quantile) / alpha_q(node_macro, quantile))
                ) * duration_period(year_all)
            )
            + SUM((year_all, quantile)$last_period(year_all),
                udf(node_macro, year_all) * (
                    LOG(CON(node_macro, year_all, quantile))
                    - beta_rc_spec_q(node_macro, quantile) * LOG(eneprice(node_macro, 'rc_spec', year_all) / 1000)
                    + beta_rc_spec_q(node_macro, quantile) * LOG(beta_rc_spec_q(node_macro, quantile) / alpha_q(node_macro, quantile))
                    - beta_rc_therm_q(node_macro, quantile) * LOG(eneprice(node_macro, 'rc_therm', year_all) / 1000)
                    + beta_rc_therm_q(node_macro, quantile) * LOG(beta_rc_therm_q(node_macro, quantile) / alpha_q(node_macro, quantile))
                    - beta_transport_q(node_macro, quantile) * LOG(eneprice(node_macro, 'transport', year_all) / 1000)
                    + beta_transport_q(node_macro, quantile) * LOG(beta_transport_q(node_macro, quantile) / alpha_q(node_macro, quantile))
                ) * duration_period(year_all)
                + 1 / finite_time_corr(node_macro, year_all)
            )
        )
    ) ;

* FOC für Haushaltskonsum pro quantile
FOC_CON(node_macro, year_all, quantile)$(macro_horizon(year_all) AND NOT macro_base_period(year_all))..
    hh_scale * udf(node_macro, year_all) * duration_period(year_all)
        / CON(node_macro, year_all, quantile)
    - SUM(year2$seq_period(year_all, year2),
        MU_CAPITAL(node_macro, year2, quantile) * duration_period(year2) * ((alpha_q(node_macro, quantile) + beta_rc_spec_q(node_macro, quantile) + beta_rc_therm_q(node_macro, quantile) + beta_transport_q(node_macro, quantile)) / alpha_q(node_macro, quantile))
    )
    - NU_TERM(node_macro, year_all, quantile)$last_period(year_all) * ((alpha_q(node_macro, quantile) + beta_rc_spec_q(node_macro, quantile) + beta_rc_therm_q(node_macro, quantile) + beta_transport_q(node_macro, quantile)) / alpha_q(node_macro, quantile))
    =G= 0 ;

* FOC für Haushaltskapital pro quantile
FOC_KAP(node_macro, year_all, quantile)$(macro_horizon(year_all) AND NOT macro_base_period(year_all))..
    MU_CAPITAL(node_macro, year_all, quantile)
    - SUM(year2$seq_period(year_all, year2),
        MU_CAPITAL(node_macro, year2, quantile)
            * ((1 - depr(node_macro)) ** duration_period(year2)
            + duration_period(year2) * (((1 + INTEREST(node_macro, year_all)) ** duration_period(year2) - 1) / duration_period(year2)))
    )
    - NU_TERM(node_macro, year_all, quantile)$last_period(year_all)
        * ((((1 + INTEREST(node_macro, year_all)) ** duration_period(year_all) - 1) / duration_period(year_all))
            - KGROW(node_macro, year_all) - depr(node_macro))
    =G= 0 ;


* Arbeitsmarktgleichgewicht: Summe Haushaltsarbeitsangebot = Aggregat
LABOR_MARKET(node_macro, year_all)$(NOT macro_base_period(year_all))..
    labor(node_macro, year_all) =E= SUM(quantile, LAB(node_macro, year_all, quantile));

* Terminalbedingung für einzelne Haushalte (letztes Jahr)
TERMINAL_HH_MCP(node_macro, year_all, quantile)$(last_period(year_all))..
    KAP(node_macro, year_all, quantile)
        * ((((1 + INTEREST(node_macro, year_all)) ** duration_period(year_all) - 1) / duration_period(year_all))
            - KGROW(node_macro, year_all) - depr(node_macro))
    + LAB(node_macro, year_all, quantile) * WAGE(node_macro, year_all)
    - eneprice(node_macro, 'rc_spec', year_all)/1000 * quantile_share(quantile) * EMIN(node_macro)
    - eneprice(node_macro, 'rc_therm', year_all)/1000 * quantile_share(quantile) * EMIN(node_macro)
    - eneprice(node_macro, 'transport', year_all)/1000 * quantile_share(quantile) * EMIN(node_macro)
    - ((alpha_q(node_macro, quantile) + beta_rc_spec_q(node_macro, quantile) + beta_rc_therm_q(node_macro, quantile) + beta_transport_q(node_macro, quantile))/alpha_q(node_macro, quantile)) * CON(node_macro, year_all, quantile)
    =E= 0 ;

* Kapitalwachstum (KGROW)
CAPITAL_GROWTH_MCP(node_macro, year_all)$(NOT macro_base_period(year_all))..
    KGROW(node_macro, year_all) =E=
        SUM(year2$(seq_period(year2, year_all)), (K(node_macro, year_all) - K(node_macro, year2)) / K(node_macro, year2));
    
* --- MCP-EMP-Block (schematisch, anpassen nach Bedarf) ---
$EMP
VI
    KAP_DYN.MU_KAP
    INV_ACC.MU_INV
    CAP_USE.LAMBDA_CAP
    PROD_FUNC.MU_PROD
    FOC_KAP_PROD.MU_MPK
    FOC_LAB_PROD.MU_MPL
    * weitere Zuordnungen
    ENERGY_ACCOUNTING_MCP.mu_ENERGY_ACCOUNTING
    ENERGY_ACCOUNTING2_MCP.mu_ENERGY_ACCOUNTING2
    ENERGY_SUPPLY_MCP.mu_ENERGY_SUPPLY
    COST_ENERGY_MCP.mu_COST_ENERGY
    AGG_KAP.MU_AGG_KAP
    AGG_CON.MU_AGG_CON
    AGG_INV.MU_AGG_INV
    CAP_USE.MU_CAP_USE
    TERMINAL_CONDITION_MCP.MU_TERM_COND
    HH_UTILITY_REP_DEF.MU_HH_UTILITY_REP
    FOC_CON.MU_FOC_CON
    FOC_KAP.MU_FOC_KAP
    LABOR_MARKET.MU_LABOR_MARKET
    TERMINAL_HH_MCP.MU_TERMINAL_HH
    CAPITAL_GROWTH_MCP.MU_KGROW
/

* Modell- und Solve-Block folgt nach vollständiger Definition aller Blöcke
* ...
