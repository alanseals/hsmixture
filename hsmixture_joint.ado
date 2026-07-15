*! version 2.4.0  14jul2026
*! Heckman-Singer Joint Timing-of-Events Model
*! Two-equation discrete-time hazard with correlated unobserved heterogeneity.
*! Supports two factor structures via factor() option:
*!   common   - one shared loading lambda (manuscript-style one-factor MPH)
*!   separate - free loadings lambda_T, lambda_Y (default; admits opposite
*!              signs)
*!
*! Authors: Jonghoon Park and R. Alan Seals (Auburn University)
*!
*! References:
*!   Abbring, J.H. and G.J. van den Berg. 2003. The nonparametric
*!     identification of treatment effects in duration models.
*!     Econometrica 71(5): 1491-1517.
*!   Heckman, J.J. and B. Singer. 1984. Econometrica 52(2): 271-320.

program hsmixture_joint, eclass sortpreserve
    version 14

    if replay() {
        if "`e(cmd)'" != "hsmixture_joint" {
            error 301
        }
        Display `0'
        exit
    }

    * Parse the compound syntax
    * Syntax: (treat_dep = treat_indep) (outcome_dep = outcome_indep, treat(varname))
    *         , id() k() [options]

    gettoken eq1 0 : 0, parse("(") match(paren)
    gettoken eq2 0 : 0, parse("(") match(paren)

    * Parse equation 1 (treatment process)
    gettoken treat_dep eq1_rest : eq1, parse("=")
    gettoken eq_sign treat_indep : eq1_rest, parse("=")
    local treat_dep = strtrim("`treat_dep'")
    local treat_indep = strtrim("`treat_indep'")

    * Parse equation 2 (outcome process with treatment indicator)
    local treat_var ""

    * Check if treat() is specified
    if regexm("`eq2'", "treat\(([^\)]+)\)") {
        local treat_var = regexs(1)
        * Remove treat() from eq2
        local eq2 = regexr("`eq2'", ",[ ]*treat\([^\)]+\)", "")
    }

    gettoken outcome_dep eq2_rest : eq2, parse("=")
    gettoken eq_sign outcome_indep : eq2_rest, parse("=")
    local outcome_dep = strtrim("`outcome_dep'")
    local outcome_indep = strtrim("`outcome_indep'")

    * Parse remaining options. PRisk() is a deprecated alias for RISKset()
    * retained for back-compat with v2.0.0 caller scripts.
    *
    * FACTor() selects the heterogeneity structure:
    *   separate (default): two free signed loadings lambda_T, lambda_Y
    *     for the treatment and outcome equations. Per-type shifts
    *     (lambda_T*v_k, lambda_Y*v_k) lie on a 1-D locus through the
    *     origin but the locus direction is free. This admits opposite
    *     signs (negative correlation between treatment-prone and
    *     outcome-prone latent types).
    *   common: one free signed loading lambda shared by both equations,
    *     so the type-k shift is (lambda*v_k, lambda*v_k). This is the
    *     classical Heckman-Singer one-factor mixed-proportional-hazards
    *     model: the same exp(eta_i) raises both hazards. Correlation
    *     between latent treatment and outcome propensity is mechanically
    *     positive.
    *
    * NSTarts(integer 7) selects the number of starting-value configurations
    * for multistart (best log-likelihood retained). Mixture surfaces are
    * multimodal, so single-start is unsafe at any K. Under factor(common)
    * the configuration grid is smaller (one fewer loading parameter);
    * the default 7 is capped to 6 in that mode.
    *
    * DIFFicult, TRace, GRADient, HESSian, TECHnique(), TOLerance(),
    * LTOLerance(), NRTOLerance() are legacy `ml model d0` options from
    * v2.0.0. They are accepted but ignored by the Mata optimize()
    * implementation, which uses BFGS with vtol=1e-10 + nrtol=1e-5; a
    * runtime note (below) makes the no-op discoverable.
    syntax [if] [in], ///
        ID(varname) ///
        [K(integer 2) ///
         FACTor(string) ///
         RISKset(varname) ///
         PRisk(varname) ///
         TIME(varname) ///
         FROM(name) ///
         ITERate(integer 100) ///
         NSTarts(integer 7) ///
         noLOG ///
         Level(cilevel) ///
         DIFFicult ///
         TRace ///
         GRADient ///
         HESSian ///
         TECHnique(string) ///
         TOLerance(real 1e-6) ///
         LTOLerance(real 1e-7) ///
         NRTOLerance(real 1e-5)]

    * Legacy ml-model-d0 optimizer options (difficult/trace/gradient/hessian/
    * technique()/tolerance()/ltolerance()/nrtolerance()) are accepted for
    * back-compatibility with v2.0.0 callers but have NO effect: estimation uses
    * Mata optimize() BFGS with fixed tolerances. Warn if a caller supplied one
    * of the detectable toggles or technique() so the silent no-op is
    * discoverable (mirrors the prisk() deprecation notice).
    if "`difficult'`trace'`gradient'`hessian'`technique'" != "" ///
        | `tolerance' != 1e-6 | `ltolerance' != 1e-7 | `nrtolerance' != 1e-5 {
        display as txt "note: legacy optimizer options (difficult, trace, gradient,"
        display as txt "  hessian, technique, tolerance, ltolerance, nrtolerance) are"
        display as txt "  accepted for back-compatibility with v2.0.0 callers and"
        display as txt "  ignored; estimation uses Mata optimize() BFGS with fixed"
        display as txt "  tolerances."
    }

    * Resolve factor() mode. Default separate for back-compat with v2.x callers.
    if "`factor'" == "" {
        local factor "separate"
    }
    if !inlist("`factor'", "common", "separate") {
        display as error ///
            "factor() must be 'common' or 'separate' (got '`factor'')"
        exit 198
    }

    * Resolve prisk/riskset alias
    if "`prisk'" != "" & "`riskset'" != "" {
        display as error "cannot specify both prisk() and riskset(); use riskset()"
        exit 198
    }
    if "`prisk'" != "" & "`riskset'" == "" {
        local riskset "`prisk'"
        display as txt "(prisk() is a deprecated alias for riskset(); please update callers)"
    }

    * Mark sample. Include riskset so that observations with missing
    * riskset values are excluded from the estimation sample rather than
    * silently entering Mata as missings (which would multiply the
    * outcome log-likelihood contribution by missing).
    marksample touse
    markout `touse' `id' `treat_dep' `outcome_dep' `treat_var' `riskset'

    * Validate inputs
    if "`treat_var'" == "" {
        display as error "treat() required in second equation"
        exit 198
    }

    if `k' < 2 {
        display as error "k() must be at least 2"
        exit 198
    }

    * Expand wildcards and exclude obs with missing covariates from touse
    * BEFORE the data-contract assertions. If a row has missing covariates it
    * will not enter the estimation sample, so it must not be checked by the
    * 0/1 / single-event / riskset assertions either.
    unab treat_indep_exp : `treat_indep'
    unab outcome_indep_exp : `outcome_indep'
    markout `touse' `treat_indep_exp' `outcome_indep_exp'

    * Within-person ordering for the row-level contract checks below (one-time
    * treatment, absorbing outcome, treated-indicator consistency). When
    * time() is supplied, order by it (validated numeric, nonmissing, unique
    * within id). Otherwise use the incoming row order, which the data
    * contract (sorted by id and time) makes meaningful. Captured HERE,
    * before any bysort below re-sorts the data: sortpreserve restores the
    * user's order on exit, but bysort changes it during execution and an
    * unstable sort can scramble within-id order.
    tempvar __hsj_ord
    if "`time'" != "" {
        capture confirm numeric variable `time'
        if _rc {
            display as error "time() variable `time' must be numeric"
            exit 198
        }
        capture assert !missing(`time') if `touse'
        if _rc {
            display as error ///
                "time() variable `time' has missing values in the estimation sample"
            exit 198
        }
        tempvar __hsj_tdup
        quietly bysort `id' `time': egen double `__hsj_tdup' = total(`touse')
        capture assert `__hsj_tdup' <= 1 if `touse'
        if _rc {
            display as error ///
                "time() variable `time' has duplicate values within id in the estimation sample"
            display as error ///
                "  each person-period must have a distinct time value"
            exit 198
        }
        drop `__hsj_tdup'
        quietly gen double `__hsj_ord' = `time'
    }
    else {
        * No time() supplied: the row contracts below read within-person
        * order from physical row order, which is what the documented data
        * contract (sorted by id and time) delivers. Time order cannot be
        * verified without a time variable, but a dataset whose sort key
        * does not even begin with `id' has almost certainly been reordered
        * (merge, append, expand), and then row order carries no time
        * information -- which would let the checks below both miss real
        * violations and flag phantom ones. Warn in that case. This is a
        * note rather than an error because the physical order may still be
        * correct and the command cannot prove otherwise.
        local __hsj_sortkey : sortedby
        local __hsj_s1 = word("`__hsj_sortkey'", 1)
        if "`__hsj_s1'" != "`id'" {
            display as txt ///
                "note: the data are not sorted by `id'. The data-contract checks read"
            display as txt ///
                "  within-person order from the current row order. If rows are not in time"
            display as txt ///
                "  order within `id', sort the data or supply time(varname); otherwise the"
            display as txt ///
                "  checks may miss a malformed panel. The likelihood itself is unaffected."
        }
        quietly gen double `__hsj_ord' = _n
    }

    * Stage 3: data-contract validation. Errors (not warnings) on violations
    * because these break the likelihood. Risk-set rows are dropped via
    * markout above; here we check the values that survive.
    capture assert inlist(`treat_dep', 0, 1) if `touse'
    if _rc {
        display as error ///
            "treatment event variable `treat_dep' must be 0/1 in the estimation sample"
        exit 198
    }
    capture assert inlist(`outcome_dep', 0, 1) if `touse'
    if _rc {
        display as error ///
            "outcome event variable `outcome_dep' must be 0/1 in the estimation sample"
        exit 198
    }
    capture assert inlist(`treat_var', 0, 1) if `touse'
    if _rc {
        display as error ///
            "treatment indicator `treat_var' must be 0/1 in the estimation sample"
        exit 198
    }
    if "`riskset'" != "" {
        capture assert inlist(`riskset', 0, 1) if `touse'
        if _rc {
            display as error "riskset() variable `riskset' must be 0/1 in the estimation sample"
            exit 198
        }
    }

    * Stage 3 (extended): duration-data structural validation. The likelihood
    * assumes one-time treatment events (each person becomes treated at most
    * once) and an absorbing outcome (each person experiences the outcome at
    * most once). Violations of these assumptions will distort delta and the
    * mixture parameters in subtle ways. Errors here protect users from
    * silently feeding bad person-period panels into the estimator.
    tempvar __hsj_evcount
    quietly bysort `id': egen `__hsj_evcount' = total(`treat_dep') if `touse'
    capture assert `__hsj_evcount' <= 1 if `touse'
    if _rc {
        display as error ///
            "treatment event variable `treat_dep' has more than one event for some id"
        display as error ///
            "  hsmixture_joint requires one-time treatment timing (each person treated at most once)"
        exit 198
    }
    drop `__hsj_evcount'

    tempvar __hsj_evcount
    quietly bysort `id': egen `__hsj_evcount' = total(`outcome_dep') if `touse'
    capture assert `__hsj_evcount' <= 1 if `touse'
    if _rc {
        display as error ///
            "outcome event variable `outcome_dep' has more than one event for some id"
        display as error ///
            "  hsmixture_joint requires an absorbing outcome (each person experiences the outcome at most once)"
        exit 198
    }
    drop `__hsj_evcount'

    * If riskset is specified, require treat_event=1 only when riskset=1.
    if "`riskset'" != "" {
        capture assert `riskset' == 1 if `treat_dep' == 1 & `touse'
        if _rc {
            display as error ///
                "treatment events occur outside the treatment risk set (`treat_dep'=1 with `riskset'=0)"
            display as error ///
                "  treatment can only happen when the person is at risk (`riskset'=1)"
            exit 198
        }
    }

    * ----------------------------------------------------------------
    * Stage 3 (row contracts, v2.4.0). The one-time-treatment and
    * absorbing-outcome assumptions are contracts on ROWS, not just on
    * event counts: the likelihood sums a contribution over every
    * surviving estimation row, so rows that should not be at risk must
    * not be present (outcome) or must be excluded via riskset (treatment).
    * These checks enforce that, using `__hsj_ord' for within-person order.
    * ----------------------------------------------------------------
    tempvar __hsj_evT __hsj_evY
    quietly bysort `id': egen double `__hsj_evT' = ///
        min(cond(`treat_dep' == 1 & `touse', `__hsj_ord', .))
    quietly bysort `id': egen double `__hsj_evY' = ///
        min(cond(`outcome_dep' == 1 & `touse', `__hsj_ord', .))

    * (a) Absorbing outcome: no estimation rows after the outcome event.
    capture assert missing(`__hsj_evY') | `__hsj_ord' <= `__hsj_evY' if `touse'
    if _rc {
        display as error ///
            "estimation-sample rows occur after the outcome event for some id"
        display as error ///
            "  hsmixture_joint assumes an absorbing outcome: a person's rows must stop"
        display as error ///
            "  at the outcome event row. Rows after the event enter the likelihood as"
        display as error ///
            "  spurious at-risk periods and bias delta and the mixture parameters."
        display as error ///
            "  Drop post-event rows before estimating."
        if "`time'" != "" {
            display as error "  Ordering used: time(`time')."
        }
        else {
            display as error ///
                "  Ordering used: current row order within id. If the data are not sorted"
            display as error ///
                "  by id and time, sort them first or supply the time() option."
        }
        exit 198
    }

    * (b) One-time treatment: person-periods after the treatment event must
    * not re-enter the treatment-equation likelihood as at-risk periods.
    if "`riskset'" != "" {
        capture assert `riskset' == 0 ///
            if `touse' & !missing(`__hsj_evT') & `__hsj_ord' > `__hsj_evT'
        if _rc {
            display as error ///
                "riskset() rows occur after the treatment event (`riskset'=1 on post-event rows)"
            display as error ///
                "  hsmixture_joint models one-time treatment timing: person-periods after"
            display as error ///
                "  the treatment event are no longer at risk and must have `riskset'=0."
            exit 198
        }
    }
    else {
        * Without riskset(), the treatment equation treats EVERY estimation
        * row as at risk. That is only correct when no treated person has
        * rows after the treatment event; otherwise post-event rows add
        * spurious survival terms (the classic silent mis-specification).
        capture assert missing(`__hsj_evT') | `__hsj_ord' <= `__hsj_evT' if `touse'
        if _rc {
            display as error "riskset() is required for this data"
            display as error ///
                "  Some persons have estimation rows after their treatment event. Without"
            display as error ///
                "  a treatment risk set those rows enter the treatment-equation likelihood"
            display as error ///
                "  as spurious at-risk periods and bias the loadings and mixture. Supply a"
            display as error ///
                "  0/1 indicator equal to 1 up to and including the treatment event, e.g."
            display as error ///
                "      gen byte trisk = (`treat_var' == 0)"
            display as error ///
                "  (valid when `treat_var' switches on the period after the event) and pass"
            display as error ///
                "  riskset(trisk)."
            exit 198
        }
    }

    * (c) Treatment-indicator consistency: treated is an absorbing state and
    * must not lead its own event. Two conventions are accepted -- treated
    * switches on in the event period (same-period) or in the following
    * period (lagged, the no-anticipation convention) -- but it must never
    * be 1 before the event row and must never revert to 0.
    tempvar __hsj_t1min __hsj_t0max
    quietly bysort `id': egen double `__hsj_t1min' = ///
        min(cond(`treat_var' == 1 & `touse', `__hsj_ord', .))
    quietly bysort `id': egen double `__hsj_t0max' = ///
        max(cond(`treat_var' == 0 & `touse', `__hsj_ord', .))

    capture assert missing(`__hsj_t1min') | missing(`__hsj_t0max') | ///
        `__hsj_t0max' < `__hsj_t1min' if `touse'
    if _rc {
        display as error ///
            "treatment indicator `treat_var' reverts from 1 to 0 within id"
        display as error ///
            "  hsmixture_joint models an absorbing treatment state: once `treat_var'"
        display as error ///
            "  switches to 1 it must remain 1 for the person's remaining rows."
        exit 198
    }

    capture assert missing(`__hsj_evT') | missing(`__hsj_t1min') | ///
        `__hsj_t1min' >= `__hsj_evT' if `touse'
    if _rc {
        display as error ///
            "treatment indicator `treat_var' equals 1 before the treatment event for some id"
        display as error ///
            "  `treat_var' may switch on in the event period or later (same-period and"
        display as error ///
            "  lagged conventions are both accepted), never before the event."
        exit 198
    }

    * Persons who switch to treated with no treatment event row in the
    * estimation sample: their treatment timing is not modeled (the
    * treatment equation sees them as censored). Warn rather than error --
    * an intentionally excluded event row (e.g., an event outside an
    * eligibility window) is a legitimate design, an accidentally dropped
    * one is not. Always-treated persons with no event row (left-truncated
    * entry) are accepted silently; see the help file.
    tempvar __hsj_tag
    quietly egen byte `__hsj_tag' = tag(`id') if `touse'
    quietly count if `__hsj_tag' == 1 & !missing(`__hsj_t1min') & ///
        !missing(`__hsj_t0max') & missing(`__hsj_evT')
    if r(N) > 0 {
        display as txt "note: " as res r(N) as txt ///
            " person(s) switch to `treat_var'=1 with no `treat_dep'=1 row in the"
        display as txt ///
            "  estimation sample. Their treatment timing does not contribute to the"
        display as txt ///
            "  treatment equation (censored there); verify the event rows were"
        display as txt "  excluded intentionally."
    }
    drop `__hsj_ord' `__hsj_evT' `__hsj_evY' `__hsj_t1min' `__hsj_t0max' `__hsj_tag'

    * Compile Mata functions (silently skips if already loaded)
    capture findfile hsmixture_joint_mata.do
    if _rc {
        display as error "hsmixture_joint_mata.do not found on adopath"
        exit 601
    }
    capture quietly do "`r(fn)'"

    * ----------------------------------------------------------------
    * Set globals for likelihood evaluator (cleaned up on all exit paths)
    * ----------------------------------------------------------------
    global HSJ_K = `k'
    global HSJ_id "`id'"
    global HSJ_treat_event "`treat_dep'"
    global HSJ_outcome_event "`outcome_dep'"
    global HSJ_treat "`treat_var'"
    global HSJ_riskset "`riskset'"
    global HSJ_factor "`factor'"
    if "`log'" == "nolog" {
        global HSJ_nolog "1"
    }
    else {
        global HSJ_nolog ""
    }

    * Wrap estimation in capture for guaranteed cleanup on break/error
    capture noisily {

    * ----------------------------------------------------------------
    * Fit initial GLMs to get starting values and expand factor vars
    * ----------------------------------------------------------------
    if "`log'" != "nolog" {
        display _n as txt "Obtaining starting values from separate GLM models..."
    }

    * Treatment equation
    if "`riskset'" != "" {
        quietly glm `treat_dep' `treat_indep' if `touse' & `riskset' == 1, ///
            family(binomial) link(cloglog) nolog
    }
    else {
        quietly glm `treat_dep' `treat_indep' if `touse', ///
            family(binomial) link(cloglog) nolog
    }
    tempname b_treat
    matrix `b_treat' = e(b)

    * `treat_indep_exp' and `outcome_indep_exp' were expanded above (before
    * the data-contract asserts) so touse already excludes obs with missing
    * covariates.
    local treat_vars_mata "`treat_indep_exp'"

    * Outcome equation
    quietly glm `outcome_dep' `outcome_indep' if `touse', ///
        family(binomial) link(cloglog) nolog
    tempname b_outcome
    matrix `b_outcome' = e(b)

    local outcome_vars_mata "`outcome_indep_exp'"

    * Store variable names in globals for Mata
    global HSJ_treat_vars_mata "`treat_vars_mata'"
    global HSJ_outcome_vars_mata "`outcome_vars_mata'"

    * Get starting value for delta from simple cloglog with treatment
    quietly glm `outcome_dep' `treat_var' `outcome_indep' if `touse', ///
        family(binomial) link(cloglog) nolog
    local delta_init = _b[`treat_var']

    * ----------------------------------------------------------------
    * Build multistart starting-value configurations
    * ----------------------------------------------------------------
    * Layout under v_2=1 normalization:
    *   factor(separate): [beta_T | beta_Y | delta | lambda_T | lambda_Y |
    *                      v_3..v_K | eta_2..eta_K]
    *   factor(common):   [beta_T | beta_Y | delta | lambda |
    *                      v_3..v_K | eta_2..eta_K]
    * Loadings are signed real-line parameters. v_3..v_K (when K >= 3) are
    * free mass points. eta_k feeds the softmax mixture.
    *
    * Multistart probes regions of the loading space:
    *   separate: 7 configs. Configs 1-6 vary (lambda_T, lambda_Y) magnitudes
    *     and (for K>=3) v_3, with balanced initial mixture. Config 7 probes
    *     the small-loading + asymmetric-mixture region.
    *   common: 6 configs. The loading is one-dimensional, so we don't need
    *     the asymmetric (lambda_T != lambda_Y) probes from configs 3-4.
    *
    * Best log-likelihood across all starts wins.

    if "`factor'" == "common" {
        * Common-loading configs. lambda is the single shared loading.
        local sv_l_1   = 0.3
        local sv_v3_1  = 2.0
        local sv_eta2_1 = 0

        local sv_l_2   = 0.5
        local sv_v3_2  = 2.0
        local sv_eta2_2 = 0

        local sv_l_3   = 1.0
        local sv_v3_3  = 2.5
        local sv_eta2_3 = 0

        local sv_l_4   = 1.0
        local sv_v3_4  = 3.0
        local sv_eta2_4 = 0

        local sv_l_5   = 1.5
        local sv_v3_5  = 4.0
        local sv_eta2_5 = 0

        * Config 6: small-loading, asymmetric-mixture probe
        local sv_l_6   = 0.05
        local sv_v3_6  = 2.0
        local sv_eta2_6 = 1.9

        local n_configs = min(`nstarts', 6)
        if `n_configs' < 1 local n_configs = 1

        forvalues s = 1/`n_configs' {
            tempname b0_`s'
            matrix `b0_`s'' = `b_treat', `b_outcome', `delta_init'

            * Single shared loading
            matrix `b0_`s'' = `b0_`s'', `sv_l_`s''

            * v_3..v_K (empty for K=2)
            forvalues j = 3/`k' {
                local v_init = `sv_v3_`s'' + (`j' - 3)
                matrix `b0_`s'' = `b0_`s'', `v_init'
            }

            * eta_2..eta_K
            matrix `b0_`s'' = `b0_`s'', `sv_eta2_`s''
            forvalues j = 3/`k' {
                matrix `b0_`s'' = `b0_`s'', 0
            }
        }
    }
    else {
        * Separate-loading configs (current default). Per-start
        * (lambda_T, lambda_Y, v_3, eta_2) values; v_4..v_K built off v_3.
        local sv_lT_1   = 0.3
        local sv_lY_1   = 0.3
        local sv_v3_1   = 2.0
        local sv_eta2_1 = 0

        local sv_lT_2   = 0.5
        local sv_lY_2   = 0.5
        local sv_v3_2   = 2.0
        local sv_eta2_2 = 0

        local sv_lT_3   = 1.0
        local sv_lY_3   = 0.5
        local sv_v3_3   = 2.0
        local sv_eta2_3 = 0

        local sv_lT_4   = 0.5
        local sv_lY_4   = 1.0
        local sv_v3_4   = 2.5
        local sv_eta2_4 = 0

        local sv_lT_5   = 1.0
        local sv_lY_5   = 1.0
        local sv_v3_5   = 3.0
        local sv_eta2_5 = 0

        local sv_lT_6   = 1.5
        local sv_lY_6   = 1.5
        local sv_v3_6   = 4.0
        local sv_eta2_6 = 0

        * Config 7: small-loading, asymmetric-mixture probe
        local sv_lT_7   = 0.04
        local sv_lY_7   = 0.06
        local sv_v3_7   = 2.0
        local sv_eta2_7 = 1.9

        local n_configs = min(`nstarts', 7)
        if `n_configs' < 1 local n_configs = 1

        forvalues s = 1/`n_configs' {
            tempname b0_`s'
            matrix `b0_`s'' = `b_treat', `b_outcome', `delta_init'

            * lambda_T, lambda_Y (signed real, no exp transform)
            matrix `b0_`s'' = `b0_`s'', `sv_lT_`s'', `sv_lY_`s''

            * v_3..v_K (empty for K=2)
            forvalues j = 3/`k' {
                local v_init = `sv_v3_`s'' + (`j' - 3)
                matrix `b0_`s'' = `b0_`s'', `v_init'
            }

            * eta_2..eta_K
            matrix `b0_`s'' = `b0_`s'', `sv_eta2_`s''
            forvalues j = 3/`k' {
                matrix `b0_`s'' = `b0_`s'', 0
            }
        }
    }

    * Use user-supplied starting values if provided (replaces config 1 only)
    if "`from'" != "" {
        matrix `b0_1' = `from'
    }

    * Display header
    if "`log'" != "nolog" {
        display _n as txt "Joint Timing-of-Events Model with Heckman-Singer Heterogeneity"
        display as txt "K = `k' mass points (v_1 = 0, v_2 = 1 normalized)"
        if "`factor'" == "common" {
            display as txt "Factor structure: common (one shared loading lambda)"
        }
        else {
            display as txt "Factor structure: separate (lambda_T, lambda_Y free)"
        }
        display as txt "Multistart: `n_configs' starting configurations"
        if "`riskset'" != "" {
            display as txt "Treatment risk set: `riskset' == 1"
        }
        display as txt "{hline 70}"
        display as txt "Treatment equation: `treat_dep' = `treat_indep'"
        display as txt "Outcome equation:   `outcome_dep' = `outcome_indep' + delta*`treat_var'"
        display as txt "{hline 70}"
    }

    * Ensure data is sorted by id for Mata panelsetup
    sort `id'

    * Initialize Mata cache
    capture mata: _hsj_cleanup()

    * ----------------------------------------------------------------
    * Run optimization with Mata optimize()
    * ----------------------------------------------------------------

    local best_ll = .
    local best_converged = 0
    local best_clip = .
    local best_v_scaffold = 0
    local n_finite = 0
    local n_conv_bfgs = 0
    local n_aborted = 0

    forvalues s = 1/`n_configs' {
        if "`log'" != "nolog" {
            if "`factor'" == "common" {
                if `k' >= 3 {
                    display _n as txt "--- Starting configuration `s'/`n_configs': " ///
                        "lambda=`sv_l_`s'', " ///
                        "v_3=`sv_v3_`s'', eta_2=`sv_eta2_`s'' ---"
                }
                else {
                    display _n as txt "--- Starting configuration `s'/`n_configs': " ///
                        "lambda=`sv_l_`s'', " ///
                        "eta_2=`sv_eta2_`s'' ---"
                }
            }
            else {
                if `k' >= 3 {
                    display _n as txt "--- Starting configuration `s'/`n_configs': " ///
                        "lambda_T=`sv_lT_`s'', lambda_Y=`sv_lY_`s'', " ///
                        "v_3=`sv_v3_`s'', eta_2=`sv_eta2_`s'' ---"
                }
                else {
                    display _n as txt "--- Starting configuration `s'/`n_configs': " ///
                        "lambda_T=`sv_lT_`s'', lambda_Y=`sv_lY_`s'', " ///
                        "eta_2=`sv_eta2_`s'' ---"
                }
            }
        }

        matrix __hsj_start = `b0_`s''
        scalar __hsj_has_result = 0
        scalar __hsj_converged = 0

        capture noisily mata: _hsj_run_optimize("__hsj_start", `iterate')
        local __mata_rc = _rc

        * Optimizer-abort recovery (v2.4.0). Mata's optimize() can abort
        * outright ("could not calculate numerical derivatives -- discontinuous
        * region with missing values encountered") when a numeric-derivative
        * probe degenerates near a flat optimum. Under Stata/MP the arithmetic
        * path varies run to run (parallel reduction order), so the abort is a
        * knife-edge event: the same start on the same data can complete on one
        * run and abort on the next (observed 14jul2026: identical certification
        * runs completed in 73 and 85 iterations or aborted at 62). The abort is
        * NOT a property of the likelihood, which is finite everywhere (missing
        * inputs clip; ln(0) floors at -700). Recovery: reset the Mata cache
        * (an aborted call can leave module state behind; a rebuild costs one
        * pass over the data) and retry this start once from a deterministically
        * jittered vector -- 0.1% scaling puts BFGS on a materially different
        * arithmetic path, clearing the knife edge without changing which
        * basin the start explores. A second abort falls through to the normal
        * "no usable result" handling, so one bad start can never kill the run.
        if `__mata_rc' != 0 {
            local ++n_aborted
            if "`log'" != "nolog" {
                display as txt ///
                    "  Configuration `s': optimizer aborted (rc=`__mata_rc'); retrying from a jittered start"
            }
            capture mata: _hsj_cleanup()
            matrix __hsj_start = `b0_`s'' * 1.001
            scalar __hsj_has_result = 0
            scalar __hsj_converged = 0
            capture noisily mata: _hsj_run_optimize("__hsj_start", `iterate')
        }

        * Check if Mata posted a result (finite LL)
        capture confirm scalar __hsj_has_result
        if _rc == 0 & scalar(__hsj_has_result) == 1 {
            local ll_s = scalar(__hsj_ll)
            local ++n_finite
            if scalar(__hsj_converged) == 1 {
                local ++n_conv_bfgs
            }
            if "`log'" != "nolog" {
                display as txt "  Log-likelihood: " %12.4f `ll_s'
            }

            * Selection rule: the start with the best FINITE log-likelihood
            * wins, regardless of its convergence flag. A higher-LL point is
            * the better MLE candidate even when BFGS stopped short there;
            * preferring a converged-but-lower-LL start would report a local
            * mode as if it were the optimum. The selected start's
            * convergence status is reported honestly via e(converged) /
            * e(converged_bfgs) and the strict-convergence warning.
            if `ll_s' > `best_ll' | `best_ll' == . {
                local best_ll = `ll_s'
                local best_start = `s'
                local best_converged = scalar(__hsj_converged)
                local best_ic = scalar(__hsj_ic)
                * Clip-hit count at this start's optimum (see Mata
                * _hsj_count_clips).
                local best_clip = .
                capture confirm scalar __hsj_clip
                if _rc == 0 {
                    local best_clip = scalar(__hsj_clip)
                }
                matrix __hsj_best_b = __hsj_b
                local best_v_scaffold = 0
                capture confirm matrix __hsj_V
                if _rc == 0 {
                    matrix __hsj_best_V = __hsj_V
                }
                else {
                    * Hessian inversion failed. Post a finite scaffold V so
                    * `ereturn post` succeeds and callers can still read
                    * point estimates from e(b). A missing-value V (J(np,np,.))
                    * is rejected by Stata's `ereturn post`, which would
                    * leave the entire fit unreadable. The eigenvalue check
                    * below sets e(v_pd)=0, e(v_scaffold)=1 records that this
                    * V is a placeholder, and the Display block warns
                    * prominently; callers should gate on e(converged)==1
                    * (which requires V positive definite) before trusting SEs.
                    local _np = colsof(__hsj_best_b)
                    matrix __hsj_best_V = I(`_np') * 1e-20
                    local best_v_scaffold = 1
                }
                capture confirm matrix __hsj_g
                if _rc == 0 {
                    matrix __hsj_best_g = __hsj_g
                }
                else {
                    capture matrix drop __hsj_best_g
                }
            }
        }
        else {
            if "`log'" != "nolog" {
                display as txt ///
                    "  Configuration `s' produced no usable result (non-finite log-likelihood)"
            }
        }

        capture matrix drop __hsj_start __hsj_b __hsj_V __hsj_g
        capture scalar drop __hsj_ll __hsj_has_result __hsj_converged __hsj_ic __hsj_clip
    }

    * ----------------------------------------------------------------
    * Check that at least one configuration produced results
    * ----------------------------------------------------------------

    if `best_ll' == . {
        display as error ///
            "no starting configuration produced a finite log-likelihood"
        display as error ///
            "  (this is an optimization failure, not a convergence warning; try"
        display as error ///
            "  different starting values via from() or a simpler specification)"
        exit 430
    }

    if "`log'" != "nolog" {
        display _n as txt "Best result from configuration `best_start'" ///
            " (log-lik: " %12.4f `best_ll' ")"
        display as txt "Starts: `n_finite'/`n_configs' returned a finite " ///
            "log-likelihood; `n_conv_bfgs' satisfied the BFGS criterion"
    }

    * ----------------------------------------------------------------
    * Post results to ereturn
    * ----------------------------------------------------------------

    * Equation labels use the actual dependent-variable names so that
    * post-estimation `_b[depvar:varname]` calls match what the user wrote.
    * v_2 is fixed at 1 (scale normalization) and is not a free parameter,
    * so v starts at v_3 here.
    local colnames ""
    foreach v of local treat_indep_exp {
        local colnames "`colnames' `treat_dep':`v'"
    }
    local colnames "`colnames' `treat_dep':_cons"
    foreach v of local outcome_indep_exp {
        local colnames "`colnames' `outcome_dep':`v'"
    }
    local colnames "`colnames' `outcome_dep':_cons"
    if "`factor'" == "common" {
        local colnames "`colnames' delta:_cons lambda:_cons"
    }
    else {
        local colnames "`colnames' delta:_cons lambda_T:_cons lambda_Y:_cons"
    }
    forvalues j = 3/`k' {
        local colnames "`colnames' v_`j':_cons"
    }
    forvalues j = 2/`k' {
        local colnames "`colnames' eta_`j':_cons"
    }

    matrix colnames __hsj_best_b = `colnames'
    matrix colnames __hsj_best_V = `colnames'
    matrix rownames __hsj_best_V = `colnames'

    quietly count if `touse'
    local N_obs = r(N)

    * Person-level count for AIC/BIC. The IID unit in this mixture model is
    * the person; the latent type is integrated out at the person level (see
    * Mata _hsj_compute_ll). Using person-period N would inflate the BIC
    * penalty. e(N) follows Stata's row-count convention; e(N_persons) is
    * the statistically appropriate denominator for cross-K comparisons.
    tempvar __hsj_pid_tag
    quietly egen byte `__hsj_pid_tag' = tag(`id') if `touse'
    quietly count if `__hsj_pid_tag' == 1
    local N_persons = r(N)
    drop `__hsj_pid_tag'

    ereturn post __hsj_best_b __hsj_best_V, obs(`N_obs')
    ereturn scalar ll = `best_ll'
    ereturn scalar N_persons = `N_persons'

    * Clean up Mata cache
    capture mata: _hsj_cleanup()

    * Strict-convergence diagnostics. e(converged)=1 only if (a) BFGS reported
    * convergence, (b) the relative gradient |grad|/(1+|LL|) < 1e-5, and
    * (c) the variance matrix is positive definite. The BFGS flag alone is
    * preserved as e(converged_bfgs).
    *
    * The relative metric is used because absolute |grad| does not scale with
    * sample size. For LL on the order of 10^4 the locally flat optimum has
    * |grad| ~ 5e-3 by Mata's BFGS at vtol=1e-10/nrtol=1e-5; tightening
    * further sends BFGS into a contraction-reset loop. The relative metric
    * matches the form of Mata's own nrtol stopping criterion.
    *
    * Callers gating on e(converged) will therefore distinguish
    * a clean interior optimum from a boundary "convergence" that satisfies
    * BFGS's relative-LL stopping rule but is not a stationary point.
    local strict_converged = `best_converged'
    local grad_norm = .
    local rel_grad = .
    local v_pd = 0

    capture confirm matrix __hsj_best_g
    if _rc == 0 {
        capture mata: st_numscalar("__hsj_gn", norm(st_matrix("__hsj_best_g")))
        if _rc == 0 {
            local grad_norm = scalar(__hsj_gn)
            local rel_grad = `grad_norm' / (1 + abs(`best_ll'))
            if `rel_grad' > 1e-5 local strict_converged = 0
        }
        else {
            local strict_converged = 0
        }
        capture scalar drop __hsj_gn
    }
    else {
        local strict_converged = 0
    }

    * Scale-relative positive-definiteness (v2.4.0). The old absolute floor
    * (min eigenvalue > 1e-8) made the verdict scale-dependent and, worse,
    * nondeterministic: on the opposite-signs certification the SAME fit on
    * identical data produced v_pd=0 at one run and v_pd=1 at the next
    * (14jul2026; Stata/MP last-bit variation moved an eigenvalue across the
    * fixed threshold), flunking a converged fit with honest SEs by coin flip.
    * Relative test: the smallest eigenvalue must be strictly positive and
    * exceed 1e-12 of the largest -- scale-free, and with ~100x margin over
    * eigenvalue rounding error (n*eps*maxeig ~ 1e-14*maxeig at n=46). The
    * I*1e-20 scaffold V, which the absolute floor happened to catch, is
    * excluded explicitly by the v_scaffold flag below.
    local v_mineig = .
    capture mata: st_numscalar("__hsj_meig", min(Re(eigenvalues(st_matrix("e(V)")))))
    local __rc_eig = _rc
    capture mata: st_numscalar("__hsj_xeig", max(Re(eigenvalues(st_matrix("e(V)")))))
    if `__rc_eig' == 0 & _rc == 0 {
        if !missing(scalar(__hsj_meig)) local v_mineig = scalar(__hsj_meig)
        if !missing(scalar(__hsj_meig)) & !missing(scalar(__hsj_xeig)) {
            if scalar(__hsj_meig) > 0 & ///
                scalar(__hsj_meig) > 1e-12 * scalar(__hsj_xeig) {
                local v_pd = 1
            }
        }
    }
    if `best_v_scaffold' local v_pd = 0
    if !`v_pd' local strict_converged = 0
    capture scalar drop __hsj_meig __hsj_xeig

    * Store results
    ereturn local cmd "hsmixture_joint"
    ereturn local cmdline "hsmixture_joint `0'"
    ereturn local treat_depvar "`treat_dep'"
    ereturn local outcome_depvar "`outcome_dep'"
    ereturn local treat_var "`treat_var'"
    ereturn local idvar "`id'"
    ereturn local riskset_var "`riskset'"
    ereturn local factor "`factor'"
    ereturn scalar K = `k'
    ereturn scalar k = colsof(e(b))
    * Force estimates stats / IC machinery to use the *design* parameter
    * count rather than rank(V). When the Hessian is rank-deficient
    * (e.g., a mass-point parameter pinned to a boundary, missing SE),
    * Stata's default would compute df = rank(V) < colsof(e(b)) and
    * report an AIC/BIC that disagrees with the package's own display.
    ereturn scalar rank = colsof(e(b))
    ereturn scalar df_m = colsof(e(b))
    ereturn scalar converged = `strict_converged'
    ereturn scalar converged_bfgs = `best_converged'
    ereturn scalar grad_norm = `grad_norm'
    ereturn scalar rel_grad = `rel_grad'
    ereturn scalar v_pd = `v_pd'
    ereturn scalar v_scaffold = `best_v_scaffold'
    ereturn scalar v_mineig = `v_mineig'
    ereturn scalar clip_hits = `best_clip'
    ereturn scalar ic = `best_ic'
    ereturn scalar n_starts = `n_configs'
    ereturn scalar n_finite = `n_finite'
    ereturn scalar n_bfgs_conv = `n_conv_bfgs'
    ereturn scalar n_aborted = `n_aborted'
    ereturn scalar best_start = `best_start'

    * Post gradient at best-start optimum (consumed by
    * hsmixture_joint_postestimation, convergence)
    capture confirm matrix __hsj_best_g
    if _rc == 0 {
        matrix colnames __hsj_best_g = `colnames'
        ereturn matrix gradient = __hsj_best_g
    }

    * Compute and store transformed parameters. Loadings are signed real-line
    * (Heckman-Singer factor loadings under v_2=1 scale normalization). The
    * v2.0.0 sigma_T/sigma_Y/sigma_P/sigma_D scalars are preserved as aliases
    * pointing to the same value, since under v_2=1 the "factor loading" and
    * the "type-2 risk shift" coincide.
    *
    * Under factor(common), the single shared loading lambda is exposed
    * directly and also as both _T and _Y aliases (and corresponding sigma
    * aliases). Existing v2.x callers reading e(lambda_T) and e(lambda_Y)
    * therefore work without modification — they just see equal values.
    tempname lambda_T lambda_Y delta_est hr
    scalar `delta_est' = _b[/delta]
    if "`factor'" == "common" {
        scalar `lambda_T' = _b[/lambda]
        scalar `lambda_Y' = `lambda_T'
        ereturn scalar lambda = `lambda_T'
    }
    else {
        scalar `lambda_T' = _b[/lambda_T]
        scalar `lambda_Y' = _b[/lambda_Y]
    }
    scalar `hr' = exp(`delta_est')

    ereturn scalar delta = `delta_est'
    ereturn scalar lambda_T = `lambda_T'
    ereturn scalar lambda_Y = `lambda_Y'
    ereturn scalar sigma_T = `lambda_T'
    ereturn scalar sigma_Y = `lambda_Y'
    ereturn scalar sigma_P = `lambda_T'
    ereturn scalar sigma_D = `lambda_Y'
    ereturn scalar hr = `hr'

    * Compute mixture probabilities via max-shifted softmax. Subtracting the
    * largest logit before exponentiating prevents overflow when a mixture
    * logit is large; the shift cancels in the ratio (eta_1 = 0 is the
    * reference logit).
    tempname max_eta sum_exp_eta pi
    scalar `max_eta' = 0
    forvalues j = 2/`k' {
        scalar `max_eta' = max(`max_eta', _b[/eta_`j'])
    }
    scalar `sum_exp_eta' = exp(0 - `max_eta')
    forvalues j = 2/`k' {
        scalar `sum_exp_eta' = `sum_exp_eta' + exp(_b[/eta_`j'] - `max_eta')
    }

    matrix `pi' = J(1, `k', .)
    matrix `pi'[1, 1] = exp(0 - `max_eta') / `sum_exp_eta'
    forvalues j = 2/`k' {
        matrix `pi'[1, `j'] = exp(_b[/eta_`j'] - `max_eta') / `sum_exp_eta'
    }
    ereturn matrix pi = `pi'

    * Compute mass points. v_1 = 0 and v_2 = 1 are fixed by the
    * normalization; v_3..v_K are read from the coefficient vector.
    tempname v
    matrix `v' = J(1, `k', .)
    matrix `v'[1, 1] = 0
    if `k' >= 2 {
        matrix `v'[1, 2] = 1
    }
    forvalues j = 3/`k' {
        matrix `v'[1, `j'] = _b[/v_`j']
    }
    ereturn matrix v = `v'

    * Compute CI for hazard ratio. On a fit that did not strictly converge, e(V)
    * is the I*1e-20 scaffold, so hr_lo/hr_hi collapse to a fabricated near-zero-
    * width interval. Post missing in that case so a caller reading e() without
    * checking e(converged) cannot mistake the placeholder for real precision.
    local se_delta = _se[/delta]
    local z = invnormal(1 - (1 - `level'/100)/2)
    local hr_lo = exp(`delta_est' - `z' * `se_delta')
    local hr_hi = exp(`delta_est' + `z' * `se_delta')
    if `strict_converged' {
        ereturn scalar hr_ci_lo = `hr_lo'
        ereturn scalar hr_ci_hi = `hr_hi'
    }
    else {
        ereturn scalar hr_ci_lo = .
        ereturn scalar hr_ci_hi = .
    }
    ereturn scalar level = `level'

    * Display results
    Display, level(`level')

    } // end capture noisily

    local rc = _rc

    * Guaranteed cleanup on all exit paths (normal, break, error)
    capture mata: _hsj_cleanup()
    capture macro drop HSJ_K HSJ_id HSJ_treat_event HSJ_outcome_event HSJ_treat
    capture macro drop HSJ_riskset HSJ_factor HSJ_nolog HSJ_treat_vars_mata HSJ_outcome_vars_mata
    capture matrix drop __hsj_start __hsj_b __hsj_V __hsj_g __hsj_best_b __hsj_best_V __hsj_best_g
    capture scalar drop __hsj_ll __hsj_has_result __hsj_converged __hsj_ic __hsj_clip

    * Only propagate rc when the substantive estimation didn't complete.
    * If `e(cmd)` was set to "hsmixture_joint", `ereturn post` and the
    * identifying `ereturn local cmd` both succeeded, so the caller has
    * valid e(b), e(V), e(hr), e(delta) and friends. A non-zero rc from
    * a post-post bookkeeping step (gradient confirm, scalar drop, etc.)
    * should not be propagated as an estimation failure.
    if `rc' & "`e(cmd)'" != "hsmixture_joint" {
        exit `rc'
    }
    * Estimation succeeded — explicit clean exit so a non-zero rc lingering
    * from a `capture` cleanup above does not become the program's return rc.
    exit 0
end

program Display
    syntax [, Level(cilevel)]

    local K = e(K)
    local level_val = cond("`level'" != "", `level', e(level))

    display _n as txt "{hline 78}"
    display as txt "Joint Timing-of-Events Model with Heckman-Singer Heterogeneity"
    display as txt "{hline 78}"
    display as txt "Number of obs" _col(50) "=" _col(55) %12.0fc e(N)
    display as txt "K mass points" _col(50) "=" _col(55) %12.0f e(K)
    display as txt "Log likelihood" _col(50) "=" _col(55) %12.4f e(ll)
    if "`e(factor)'" == "common" {
        display as txt "Factor structure" _col(50) "=" _col(55) ///
            "       common (one shared loading)"
    }
    else {
        display as txt "Factor structure" _col(50) "=" _col(55) ///
            "     separate (lambda_T, lambda_Y)"
    }
    if "`e(riskset_var)'" != "" {
        display as txt "Treatment risk set" _col(50) "=" _col(55) "`e(riskset_var)' == 1"
    }
    display as txt "{hline 78}"

    * Strict-convergence warning. e(converged) = 1 only when BFGS converged AND
    * relative gradient |grad|/(1+|LL|) < 1e-5 AND V is positive definite. If
    * any of those failed, the reported point estimate may not be the MLE and
    * the SEs may not be trustworthy. Surfacing this in the main output, not
    * just postestimation.
    if e(converged) == 0 {
        display _n as err "{bf:WARNING: strict convergence not achieved.}"
        if e(rel_grad) != . & e(rel_grad) > 1e-5 {
            display as err "  Relative gradient |grad|/(1+|LL|) = " ///
                %10.2e e(rel_grad) " (threshold 1e-5; |grad| = " ///
                %10.2e e(grad_norm) ")."
        }
        if e(v_pd) == 0 {
            display as err "  Variance matrix is not positive definite."
            if e(v_scaffold) == 1 {
                display as err ///
                    "  (Hessian inversion failed; a placeholder variance matrix was posted.)"
            }
        }
        if e(converged_bfgs) == 0 {
            display as err "  BFGS did not reach its own convergence criterion."
        }
        display as err "  The reported HR and CI may be unreliable. Run "
        display as err "  hsmixture_joint_postestimation, all  for details."
        display as err "{hline 78}"
    }

    * Numerical-clip diagnostic. The likelihood truncates each linear
    * predictor to [-20, 10] for overflow safety; a nonzero count at the
    * optimum means part of the fitted surface is the clipped (approximate)
    * likelihood rather than the exact cloglog likelihood.
    if e(clip_hits) > 0 & e(clip_hits) < . {
        display _n as txt "note: " as res %12.0fc e(clip_hits) ///
            as txt " linear-predictor evaluations hit the numerical bounds"
        display as txt ///
            "  [-20, 10] at the optimum. Estimates whose fitted hazards sit at the"
        display as txt ///
            "  bounds are governed by the clipped likelihood; treat them with caution."
    }

    * Display coefficient table
    ereturn display, level(`level_val')

    * Key result: Treatment effect
    local z = invnormal(1 - (1 - `level_val'/100)/2)
    local se_delta = _se[/delta]
    local hr_lo = exp(e(delta) - `z' * `se_delta')
    local hr_hi = exp(e(delta) + `z' * `se_delta')

    display _n as txt "{hline 78}"
    display as txt "{bf:Treatment Effect (Equation 1 -> Equation 2)}"
    display as txt "{hline 78}"
    display as txt "delta (log hazard ratio)" _col(30) "=" _col(35) %9.4f e(delta) ///
        _col(50) "SE = " %7.4f `se_delta'
    display as txt "Hazard Ratio" _col(30) "=" _col(35) %9.2f e(hr)
    * Print the CI only on a strictly converged fit. Otherwise e(V) is the
    * I*1e-20 scaffold and hr_lo/hr_hi are a fabricated near-zero-width band;
    * show "not available" rather than mislead a human reader (the stored
    * e(hr_ci_lo)/e(hr_ci_hi) are missing in that case for the same reason).
    if e(converged) == 1 {
        display as txt "`level_val'% CI" _col(30) "=" _col(35) ///
            "[" %5.2f `hr_lo' ", " %5.2f `hr_hi' "]"
    }
    else {
        display as txt "`level_val'% CI" _col(30) "=" _col(35) ///
            as err "not available (fit did not strictly converge)"
    }

    * Factor loadings (signed under v_2=1 normalization)
    display _n as txt "Factor Loadings (signed; v_1=0, v_2=1 normalization):"
    if "`e(factor)'" == "common" {
        display as txt "  lambda (shared)" _col(30) "=" _col(35) %9.4f e(lambda_T)
    }
    else {
        display as txt "  lambda_T (treatment)" _col(30) "=" _col(35) %9.4f e(lambda_T)
        display as txt "  lambda_Y (outcome)" _col(30) "=" _col(35) %9.4f e(lambda_Y)
    }

    * Mixture distribution
    display _n as txt "Mixture Distribution:"
    display as txt "  Type" _col(15) "Prob (pi)" _col(30) "Mass Pt (v)" ///
        _col(50) "lambda_T*v" _col(65) "lambda_Y*v"
    display as txt "  {hline 70}"
    forvalues k = 1/`K' {
        local pi_k = e(pi)[1, `k']
        local v_k = e(v)[1, `k']
        local lT_v = e(lambda_T) * `v_k'
        local lY_v = e(lambda_Y) * `v_k'
        display as txt "  `k'" _col(15) %7.4f `pi_k' _col(30) %9.4f `v_k' ///
            _col(50) %9.4f `lT_v' _col(65) %9.4f `lY_v'
    }

    * Model fit statistics. BIC denominator is e(N_persons) because the IID
    * unit in this mixture model is the person; using person-period N would
    * inflate the penalty. e(N) is preserved for Stata convention.
    local ll = e(ll)
    local k_params = e(k)
    local N_p = e(N_persons)
    local aic = -2 * `ll' + 2 * `k_params'
    local bic = -2 * `ll' + `k_params' * ln(`N_p')

    display _n as txt "Model Fit:"
    display as txt "  AIC" _col(30) "=" _col(35) %12.2f `aic'
    display as txt "  BIC" _col(30) "=" _col(35) %12.2f `bic' ///
        "  (N = " %7.0fc `N_p' " persons)"
    display as txt "  Parameters" _col(30) "=" _col(35) %12.0f `k_params'

    display as txt "{hline 78}"
end
