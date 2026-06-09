$TITLE MACRO core MCP formulation with 10 Haushaltsdezi
$EOLCOM #

* cd /home/lukas/environments/macro_uba/message_ix/message_ix/model && gams MACRO_run.gms --in=/home/lukas/environments/macro_uba/lib/python3.10/site-packages/message_ix/model/data/MsgData_MESSAGEix_ssp2_baseline_2304_add_macro.gdx --out=/home/lukas/environments/macro_uba/message_ix/message_ix/model/test_2005_mcp.gdx > MACRO_run_smart_calibration.log 2>&1

* --- MCP Variablen und Gleichungen ---

POSITIVE VARIABLES
    KAP(node, year_all, quantile)     'Kapital pro Haushalt'
    CON(node, year_all, quantile)     'Konsum pro Haushalt'
    LAB(node, year_all, quantile)      'Labor pro Haushalt'
    YE(node, sector, year_all)        'Neue Energie (POSITIVE damit YE**rho in PROD_FUNC definiert bleibt)'
    K(node, year_all)                 'Aggregiertes Kapital'
    C(node, year_all)                 'Aggregierter Konsum'
    Y(node, year_all)                 'Produktion'
    WAGE(node, year_all)              'Lohn'
    INTEREST(node, year_all)          'Zins'
    TE(node, sector, year_all)           'Gesamter Endnutzenergieinput'
    E(node, sector, year_all)            'Direkter Energiebedarf Haushalt'
    PHYSENE(node, sector, year_all)      'Physische Endnutzenergie'
    GDP(node, year_all)
;

VARIABLES
    HH_UTILITY_REP                        'Reporting only: sum of household utilities'
    KGROW(node, year_all)                  'Kapitalwachstum'
    EC(node, year_all)                   'Systemkosten'
    I(node, year_all)                 'Aggregierte Investition'
;

FREE VARIABLES
    MU_KAP(node, year_all, quantile)  'Dual zu KAP_DYN (Kapitalakkumulation, alle Perioden)'
    MU_TERMINAL_HH(node, year_all, quantile) 'Dual zu TERMINAL_HH_MCP (I_HH(T) = KAP(T)*(KGROW+depr))'
    I_HH(node, year_all, quantile)    'Investition pro Haushalt (FREE, bestimmt durch INV_ACC oder TERMINAL)'
* Keine Dualvariable für Märkte (Kapital, Arbeit): Märkte werden auf Preisvariablen gemappt, nicht auf Duals.

* Entfernt: MU_PROD, MU_MPK, MU_MPL, MU_ENERGY_ACCOUNTING, MU_ENERGY_ACCOUNTING2, MU_ENERGY_SUPPLY, MU_AGG_CON, MU_AGG_INV, MU_CAP_USE, MU_HH_UTILITY_REP, MU_COST_ENERGY
* Erklärung: Diese Duals sind nicht nötig, weil die zugehörigen Gleichungen reine Definitionen, Aggregationen oder Reporting-Gleichungen sind und keine Optimierungsbedingungen oder Märkte darstellen.
;

* --- MCP Gleichungen ---

EQUATIONS
    FOC_KAP_PROD(node, year_all)             'FOC Kapital (MPK)'
    FOC_LAB_PROD(node, year_all)             'FOC Arbeit (MPL)'
    KAP_DYN(node, year_all, quantile)     'Kapitaldynamik Haushalt (zusammengefuehrt, wie NLP EQ_KAP)'
    AGG_KAP(node, year_all)               'Aggregiertes Kapital'
    AGG_CON(node, year_all)               'Aggregierter Konsum'
    AGG_INV(node, year_all)               'Aggregierte Investition'
    CAP_USE(node, year_all)               'Aggregierte Kapitalverwendung (Y = C + I + EC)'
    PROD_FUNC(node, year_all)                'Produktionsfunktion'
    ENERGY_ACCOUNTING_MCP(node, sector, year_all)   'Energie-Bilanzgleichung (MCP)'
    ENERGY_ACCOUNTING2_MCP(node, sector, year_all)  'Haushalts-Energiebedarf (MCP)'
    ENERGY_SUPPLY_MCP(node, sector, year_all)       'Energieangebot (MCP)'
    COST_ENERGY_MCP(node, year_all)                 'Systemkosten (MCP)'
    HH_UTILITY_REP_DEF                                    'Definition der Haushaltsnutzenfunktion (repräsentativ)'
    FOC_CON(node, year_all, quantile)              'FOC fuer Haushaltskonsum pro quantile'
    FOC_KAP(node, year_all, quantile)              'FOC fuer Haushaltskapital (Euler, NOT last_period)'
    FOC_KAP_LAST(node, year_all, quantile)         'FOC fuer Haushaltskapital im letzten Jahr (=G= paired KAP)'
    MU_TERMINAL_DEF(node, year_all, quantile)      'Definition MU_TERMINAL_HH aus MU_KAP und KGROW+depr'
    INV_ACC(node, year_all, quantile)              'Haushalts-Investitionsdefinition (NOT last_period)'
    LABOR_MARKET(node, year_all)                   'Arbeitsmarktgleichgewicht'
    TERMINAL_HH_MCP(node, year_all, quantile)      'Terminalbedingung: I_HH(T) = KAP(T)*(KGROW+depr)'
    EQ_LAB(node, year_all, quantile)               'Haushaltsarbeitsangebot nach Quantil'
    CAPITAL_GROWTH_MCP(node, year_all)             'Kapitalwachstum'
    FOC_YE(node, sector, year_all)                 'FOC YE: Grenzprodukt = Grenzkosten'
;

* Produktionsfunktion
PROD_FUNC(node_active, year)$(NOT macro_base_period(year))..
    Y(node_active, year) =E=
        (ACONST(node_active) * K(node_active, year)**(rho(node_active) * kpvs(node_active))
        * labor(node_active, year)**(rho(node_active) * (1 - kpvs(node_active)))
        + SUM(sector, BCONST(node_active, sector) * YE(node_active, sector, year)**rho(node_active))
        )**(1/rho(node_active));

* Aggregierte Kapitalverwendung (Ressourcenbilanz)
CAP_USE(node_active, year)$(NOT macro_base_period(year))..
    Y(node_active, year) =E= C(node_active, year) + I(node_active, year) + EC(node_active, year);

* --- ALTE VERSION (zweiteilig, auskommentiert als Backup) ---
* I_HH(node, year_all, quantile) 'Investition pro Haushalt (FREE)' waere in VARIABLES zu deklarieren.
* INV_ACC in EQUATIONS, INV_ACC.I_HH im MODEL-Block.
*
* KAP_DYN(node_active, year, quantile)$(NOT macro_base_period(year))..
*     KAP(node_active, year, quantile) =E=
*         SUM(year2$(seq_period(year2,year)),
*             KAP(node_active, year2, quantile) * (1 - depr(node_active))**duration_period(year)
*             + duration_period(year) * I_HH(node_active, year, quantile)
*         );
*
* INV_ACC(node_active, year, quantile)$(NOT macro_base_period(year))..
*     I_HH(node_active, year, quantile) =E=
*         SUM(year2$(seq_period(year2,year)),
*             KAP(node_active, year2, quantile) * ((1 + INTEREST(node_active, year))**duration_period(year) - 1) / duration_period(year)
*             + LAB(node_active, year, quantile) * WAGE(node_active, year)
*             - eneprice(node_active, 'rc_spec', year)/1000 * quantile_share(quantile) * EMIN(node_active)
*             - eneprice(node_active, 'rc_therm', year)/1000 * quantile_share(quantile) * EMIN(node_active)
*             - eneprice(node_active, 'transport', year)/1000 * quantile_share(quantile) * EMIN(node_active)
*             - ((alpha_q(node_active, quantile) + beta_rc_spec_q(node_active, quantile) + beta_rc_therm_q(node_active, quantile) + beta_transport_q(node_active, quantile))/alpha_q(node_active, quantile)) * CON(node_active, year, quantile)
*         );
* --- ENDE ALTE VERSION ---

* =============================================================================
* KAPITALDYNAMIK UND HAUSHALTSOPTIMIERUNG
* =============================================================================
* Struktur analog NLP: CAPITAL(t, alle t) + INVESTMENT(t, t!=T) + TERMINAL_CONDITION(T)
*
* Im MCP:
*   KAP_DYN(t)      : KAP_t = KAP_{t-1}*(1-delta)^dt + dt*I_HH_t     [=E=, Dual: MU_KAP(t)]
*   INV_ACC(t!=T)   : I_HH_t = Budget_t                               [=E=, Dual: keiner -> FREE pair]
*   TERMINAL_HH(T)  : I_HH_T = KAP_T*(KGROW_T + depr)                [=E=, Dual: MU_TERMINAL_HH(T)]
*
* Daraus folgen die FOC:
*   FOC_CON(t)       : dU/dCON - MU_KAP(t)*dt*ccf >= 0               [POSITIVE CON]
*   FOC_KAP(t<T)     : MU_KAP(t) - MU_KAP(t+1)*G_{t+1} >= 0          [POSITIVE KAP]
*   FOC_KAP_LAST(T)  : MU_KAP(T) - MU_TERMINAL_HH(T)*(KGROW+depr) >= 0 [POSITIVE KAP]
*
* KAP_DYN(t) paired mit MU_KAP(t)        [FREE -> =E= always binding]
* INV_ACC(t) paired mit I_HH(t)          [FREE -> =E= always binding]
* TERMINAL_HH_MCP(T) paired mit MU_TERMINAL_HH(T) [FREE -> =E= always binding]

* Kapitalakkumulation: KAP_t = KAP_{t-1}*(1-depr)^dt + dt*I_HH_t
KAP_DYN(node_active, year, quantile)$(NOT macro_base_period(year))..
    KAP(node_active, year, quantile) =E=
        SUM(year2$(seq_period(year2,year)),
            KAP(node_active, year2, quantile) * (1 - depr(node_active))**duration_period(year)
            + duration_period(year) * I_HH(node_active, year, quantile)
        );

* Investitionsdefinition (Budget-Gleichung): nur fuer t != T
* I_HH_t = Kapitalrendite + Arbeitseinkommen - Energiekosten - Konsumausgaben
INV_ACC(node_active, year, quantile)$(NOT macro_base_period(year) AND NOT last_period(year))..
    I_HH(node_active, year, quantile) =E=
        SUM(year2$(seq_period(year2,year)),
            KAP(node_active, year2, quantile) * ((1 + INTEREST(node_active, year))**duration_period(year) - 1) / duration_period(year)
            + LAB(node_active, year, quantile) * WAGE(node_active, year)
            - eneprice(node_active, 'rc_spec', year)/1000 * quantile_share(quantile) * EMIN(node_active)
            - eneprice(node_active, 'rc_therm', year)/1000 * quantile_share(quantile) * EMIN(node_active)
            - eneprice(node_active, 'transport', year)/1000 * quantile_share(quantile) * EMIN(node_active)
            - ((alpha_q(node_active, quantile) + beta_rc_spec_q(node_active, quantile) + beta_rc_therm_q(node_active, quantile) + beta_transport_q(node_active, quantile))
               / alpha_q(node_active, quantile)) * CON(node_active, year, quantile)
        );

* Terminalbedingung (analog NLP TERMINAL_CONDITION):
* I_HH(T) = KAP(T) * (KGROW(T) + depr)  =>  genug Investition fuer Wachstum + Abschreibung
TERMINAL_HH_MCP(node_active, year, quantile)$(last_period(year))..
    I_HH(node_active, year, quantile) =E=
        KAP(node_active, year, quantile) * (KGROW(node_active, year) + depr(node_active));


* FOC Kapital (MPK = Zins)
FOC_KAP_PROD(node_active, year)$(NOT macro_base_period(year))..
    INTEREST(node_active, year) =E=
        Y(node_active, year)**(1 - rho(node_active)) * LAKL(node_active) * kpvs(node_active)
        * K(node_active, year)**(rho(node_active)*kpvs(node_active) - 1)
        * labor(node_active, year)**(rho(node_active)*(1 - kpvs(node_active)));

* FOC Arbeit (MPL = Lohn)
FOC_LAB_PROD(node_active, year)$(NOT macro_base_period(year))..
    WAGE(node_active, year) =E=
        Y(node_active, year)**(1 - rho(node_active)) * LAKL(node_active) * (1 - kpvs(node_active))
        * K(node_active, year)**(rho(node_active)*kpvs(node_active))
        * labor(node_active, year)**(rho(node_active)*(1 - kpvs(node_active)) - 1);


* FOC_YE: Grenzprodukt von YE = eneprice
FOC_YE(node_active, sector, year)$(NOT macro_base_period(year) AND h(node_active, sector) < 1)..
    eneprice(node_active, sector, year)/1000 =E=
        Y(node_active, year)**(1 - rho(node_active))
        * BCONST(node_active, sector) * YE(node_active, sector, year)**(rho(node_active) - 1);

AGG_KAP(node_active, year)$(NOT macro_base_period(year))..
    K(node_active, year) =E= SUM(quantile, KAP(node_active, year, quantile));

AGG_CON(node_active, year)$(NOT macro_base_period(year))..
    C(node_active, year) =E= SUM(quantile, CON(node_active, year, quantile));

* AGG_INV: Aggregierte Investition = Summe der Haushalts-Investitionen.
* I_HH wird durch INV_ACC (t!=T) und TERMINAL_HH_MCP (T) bestimmt.
* CAP_USE.I bestimmt I residual aus Y-C-EC, und AGG_INV verbindet I mit SUM(I_HH).
* Paarung: AGG_INV.I im MCP-Block (I ist FREE/VARIABLE)
AGG_INV(node_active, year)$(NOT macro_base_period(year))..
    I(node_active, year) =E= SUM(quantile, I_HH(node_active, year, quantile));


* --- Energiegleichungen (MCP-Form) ---

* Energie-Bilanzgleichung
ENERGY_ACCOUNTING_MCP(node_active, sector, year)$(NOT macro_base_period(year))..
    TE(node_active, sector, year) =E= YE(node_active, sector, year) + E(node_active, sector, year);

* Haushalts-Energiebedarf (quantile-spezifisch, konsistent mit INV_ACC und FOC_CON)
ENERGY_ACCOUNTING2_MCP(node_active, sector, year)$(NOT macro_base_period(year))..
    E(node_active, sector, year) =E=
        SUM(quantile,
            quantile_share(quantile) * EMIN(node_active)
            + (beta_rc_spec_q(node_active, quantile) / alpha_q(node_active, quantile))
              * CON(node_active, year, quantile) / (eneprice(node_active, 'rc_spec', year)/1000)
        ) $ sameas(sector, 'rc_spec')
    + SUM(quantile,
            quantile_share(quantile) * EMIN(node_active)
            + (beta_rc_therm_q(node_active, quantile) / alpha_q(node_active, quantile))
              * CON(node_active, year, quantile) / (eneprice(node_active, 'rc_therm', year)/1000)
        ) $ sameas(sector, 'rc_therm')
    + SUM(quantile,
            quantile_share(quantile) * EMIN(node_active)
            + (beta_transport_q(node_active, quantile) / alpha_q(node_active, quantile))
              * CON(node_active, year, quantile) / (eneprice(node_active, 'transport', year)/1000)
        ) $ sameas(sector, 'transport')
    + 0 $ (sameas(sector, 'i_spec') OR sameas(sector, 'i_therm'));

* Energieangebot: PHYSENE >= TE * aeei_factor (=G=, paired mit PHYSENE).
* Im Optimum binding: PHYSENE = TE * aeei_factor.
* Wenn slack: PHYSENE > TE*aeei -> PHYSENE auf LB (0), was TE*aeei <= 0 impliziert -> nur binding relevant.
ENERGY_SUPPLY_MCP(node_active, sector, year)$(NOT macro_base_period(year))..
    PHYSENE(node_active, sector, year) =G= TE(node_active, sector, year) * aeei_factor(node_active, sector, year);

* Systemkosten
COST_ENERGY_MCP(node_active, year)$(NOT macro_base_period(year))..
    EC(node_active, year) =E=
        (total_cost(node_active, year)/1000
        + SUM(sector, eneprice(node_active, sector, year) * 1E-3 * (PHYSENE(node_active, sector, year) - enestart(node_active, sector, year)))
        + SUM(sector, eneprice(node_active, sector, year) * 1E-3 / enestart(node_active, sector, year)
            * (PHYSENE(node_active, sector, year) - enestart(node_active, sector, year)) * (PHYSENE(node_active, sector, year) - enestart(node_active, sector, year)))
        );

* Haushaltsnutzenfunktion (repräsentativ)
HH_UTILITY_REP_DEF..
    HH_UTILITY_REP =E=
    SUM(node_active,
        1000 * (
            SUM((year, quantile)$(macro_horizon(year) AND NOT macro_base_period(year) AND NOT last_period(year)),
                udf(node_active, year) * (
                    LOG(CON(node_active, year, quantile))
                    - beta_rc_spec_q(node_active, quantile) * LOG(eneprice(node_active, 'rc_spec', year) / 1000)
                    + beta_rc_spec_q(node_active, quantile) * LOG(beta_rc_spec_q(node_active, quantile) / alpha_q(node_active, quantile))
                    - beta_rc_therm_q(node_active, quantile) * LOG(eneprice(node_active, 'rc_therm', year) / 1000)
                    + beta_rc_therm_q(node_active, quantile) * LOG(beta_rc_therm_q(node_active, quantile) / alpha_q(node_active, quantile))
                    - beta_transport_q(node_active, quantile) * LOG(eneprice(node_active, 'transport', year) / 1000)
                    + beta_transport_q(node_active, quantile) * LOG(beta_transport_q(node_active, quantile) / alpha_q(node_active, quantile))
                ) * duration_period(year)
            )
            + SUM((year, quantile)$last_period(year),
                udf(node_active, year) * (
                    LOG(CON(node_active, year, quantile))
                    - beta_rc_spec_q(node_active, quantile) * LOG(eneprice(node_active, 'rc_spec', year) / 1000)
                    + beta_rc_spec_q(node_active, quantile) * LOG(beta_rc_spec_q(node_active, quantile) / alpha_q(node_active, quantile))
                    - beta_rc_therm_q(node_active, quantile) * LOG(eneprice(node_active, 'rc_therm', year) / 1000)
                    + beta_rc_therm_q(node_active, quantile) * LOG(beta_rc_therm_q(node_active, quantile) / alpha_q(node_active, quantile))
                    - beta_transport_q(node_active, quantile) * LOG(eneprice(node_active, 'transport', year) / 1000)
                    + beta_transport_q(node_active, quantile) * LOG(beta_transport_q(node_active, quantile) / alpha_q(node_active, quantile))
                ) * duration_period(year)
                + 1 / finite_time_corr(node_active, year)
            )
        )
    ) ;

* =============================================================================
* FOC-BEDINGUNGEN (aus Lagrangian der Haushaltsoptimierung)
* =============================================================================
*
* Lagrangian:
*   L = SUM_t [ 1000*udf_t*dt * LOG(CON_t) ]
*     - SUM_t [ MU_KAP_t * (KAP_t - KAP_{t-1}*(1-d)^dt - dt*I_HH_t) ]    <- KAP_DYN
*     - SUM_{t!=T} [ I_HH_t - Budget_t(KAP_{t-1}, r_t, w_t, L_t, CON_t) ] <- INV_ACC (implicit, I_HH FREE)
*     - MU_TERM_T * (I_HH_T - KAP_T*(KGROW_T+d))                          <- TERMINAL
*
* dL/dCON_t = 1000*udf_t*dt/CON_t - MU_KAP_t*dt*ccf_q >= 0, kompl. CON_t >= 0
* dL/dKAP_t (t<T) = -MU_KAP_t + MU_KAP_{t+1}*G_{t+1} + [term from INV_ACC(t)] >= 0
*   wobei KAP_{t-1} in INV_ACC(t) mit Koeff (r_t^dt-1)/dt, Dual INV_ACC.I_HH = MU_KAP_t*dt (aus KAP_DYN)
*   -> Euler: MU_KAP_t = MU_KAP_{t+1} * G_{t+1}
* dL/dKAP_T = -MU_KAP_T + MU_TERM_T*(KGROW_T+d) >= 0, kompl. KAP_T >= 0

* FOC Konsum (alle t != base):
* CON_t erscheint in KAP_DYN(t) via I_HH_t (durch INV_ACC), Koeff = -dt*ccf_q, Dual = MU_KAP_t
FOC_CON(node_active, year, quantile)$(NOT macro_base_period(year))..
    1000 * udf(node_active, year) * duration_period(year) / CON(node_active, year, quantile)
    - MU_KAP(node_active, year, quantile)
        * duration_period(year)
        * (alpha_q(node_active, quantile) + beta_rc_spec_q(node_active, quantile) + beta_rc_therm_q(node_active, quantile) + beta_transport_q(node_active, quantile))
        / alpha_q(node_active, quantile)
    =G= 0 ;

* FOC Kapital (t != base, t < T):
* KAP(t) erscheint in:
*   - LHS KAP_DYN(t): Koeff +1, Dual MU_KAP(t)
*   - KAP_{t-1} in INV_ACC(t+1): Koeff +(r_{t+1}^dt-1)/dt, mit Dual MU_KAP(t+1)*dt (indirekt)
*   - KAP_{t-1} in KAP_DYN(t+1): Koeff +(1-d)^{dt+1}, mit Dual MU_KAP(t+1)
* Zusammen: MU_KAP(t) = MU_KAP(t+1) * [(1-d)^dt + (1+r_{t+1})^dt - 1]
FOC_KAP(node_active, year, quantile)$(NOT macro_base_period(year) AND NOT last_period(year))..
    MU_KAP(node_active, year, quantile)
    - SUM(year2$seq_period(year, year2),
        MU_KAP(node_active, year2, quantile)
            * ( (1 - depr(node_active)) ** duration_period(year2)
                + (1 + INTEREST(node_active, year2)) ** duration_period(year2) - 1 )
    )
    =G= 0 ;

* FOC Kapital im letzten Jahr:
* KAP(T) erscheint in KAP_DYN(T) [Dual MU_KAP(T)] und TERMINAL_HH_MCP(T) [Dual MU_TERMINAL_HH(T)].
* Netto-FOC: MU_KAP(T) - MU_TERMINAL_HH(T)*(KGROW+depr) >= 0, kompl. KAP(T) >= 0
FOC_KAP_LAST(node_active, year, quantile)$last_period(year)..
    MU_KAP(node_active, year, quantile)
    - MU_TERMINAL_HH(node_active, year, quantile)
        * (KGROW(node_active, year) + depr(node_active))
    =G= 0 ;

* MU_TERMINAL_HH definiert aus der Gleichgewichtsbedingung FOC_KAP_LAST (binding im Innern):
* MU_TERMINAL_HH(T) = MU_KAP(T) / (KGROW(T) + depr)
MU_TERMINAL_DEF(node_active, year, quantile)$last_period(year)..
    MU_TERMINAL_HH(node_active, year, quantile) * (KGROW(node_active, year) + depr(node_active))
    - MU_KAP(node_active, year, quantile)
    =E= 0 ;



* Haushaltsarbeitsangebot nach Quantil
EQ_LAB(node_active, year, quantile)..
    LAB(node_active, year, quantile) =E= labor(node_active, year) * quantile_share(quantile);

* Arbeitsmarktgleichgewicht: Summe Haushaltsarbeitsangebot = Aggregat
LABOR_MARKET(node_active, year)$(NOT macro_base_period(year))..
    labor(node_active, year) =E= SUM(quantile, LAB(node_active, year, quantile));

* Kapitalwachstum (KGROW)
CAPITAL_GROWTH_MCP(node_active, year)$(NOT macro_base_period(year))..
    KGROW(node_active, year) =E=
        SUM(year2$(seq_period(year2, year)), (K(node_active, year) - K(node_active, year2)) / K(node_active, year2));
    
* --- MCP-EMP-Block (schematisch, anpassen nach Bedarf) ---

MODEL MESSAGE_MACRO /
    PROD_FUNC.Y
    CAP_USE.I
*   AGG_INV.I  -- Walras: I = SUM(I_HH) folgt implizit aus Marktclearing + KAP_DYN; CAP_USE.I erzwingt Guetermarkt
    FOC_KAP_PROD.INTEREST
    FOC_LAB_PROD.WAGE
    AGG_KAP.K
    AGG_CON.C
    ENERGY_ACCOUNTING_MCP.TE
    ENERGY_ACCOUNTING2_MCP.E
    ENERGY_SUPPLY_MCP.PHYSENE
    COST_ENERGY_MCP.EC
    EQ_LAB.LAB
    TERMINAL_HH_MCP.I_HH
    CAPITAL_GROWTH_MCP.KGROW
    FOC_CON.CON
    FOC_KAP.KAP
    FOC_KAP_LAST.KAP
    MU_TERMINAL_DEF.MU_TERMINAL_HH
    KAP_DYN.MU_KAP
    INV_ACC.I_HH
    FOC_YE.YE
/ ;
* Paarungslogik:
*   PROD_FUNC.Y              : Produktionsfunktion -> Y
*   FOC_KAP_PROD.INTEREST    : MPK = r
*   FOC_LAB_PROD.WAGE        : MPL = w
*   AGG_CON.C, AGG_INV.I, AGG_KAP.K : Aggregate aus Quantilsummen
*   ENERGY_ACCOUNTING_MCP.TE : TE = YE + E -> TE (Identitaet, bestimmt TE gegeben YE,E)
*   ENERGY_ACCOUNTING2_MCP.E : E = f(C,p) -> E
*   ENERGY_SUPPLY_MCP.PHYSENE: PHYSENE >= TE*aeei -> PHYSENE (komplementaer zu PHYSENE >= 0)
*   COST_ENERGY_MCP.EC       : EC = f(PHYSENE) -> EC
*   FOC_YE.YE                : dY/dYE = dEC/dYE (via PHYSENE=TE*aeei) -> YE
*   CAP_USE.I              : Gütermarktbilanz Y = C + I + EC -> I (residuale Investition, I erscheint in der Gleichung)


*MESSAGE_MACRO.optfile = 1;

* Märkte werden auf Preisvariablen gemappt (WAGE, INTEREST), nicht auf Duals.
* Terminalbedingungen werden als eigene Gleichungen mit Duals gemappt.

* Entfernt aus dem MODEL-Block:
* - FOC_KAP_PROD, FOC_LAB_PROD, ENERGY_ACCOUNTING_MCP, ENERGY_ACCOUNTING2_MCP, ENERGY_SUPPLY_MCP, COST_ENERGY_MCP, AGG_CON, AGG_INV, HH_UTILITY_REP_DEF
* Erklärung: AGG_INV (I = sum I_HH_q) ist nicht im MCP-Block, weil I residual durch CAP_USE bestimmt wird.
* Die Konsistenz sum(I_HH_q) = I ist eine Gleichgewichtseigenschaft, die im korrekten Gleichgewicht gilt.
