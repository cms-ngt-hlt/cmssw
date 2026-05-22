#!/usr/bin/env bash

set -u

DQM_FILE="DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root"
SCRIPT_DIR="${CMSSW_BASE}/src/Validation/RecoTau/scripts"

MAKE_COMPARISON="${SCRIPT_DIR}/makeComparisonPlots.py"
MAKE_TAU_VALIDATION="${SCRIPT_DIR}/makeTauValidationPlots.py"

BASE_DIR="DQMData/Run 1/HLT/Run summary/Tau/TauValidation_DeltaR"

OUTDIR_COMPARISON="TauValidationPlots/Comparisons"
OUTDIR_SUMMARY="TauValidationPlots/Summary"
OUTDIR_RESPONSE="TauValidationPlots/Response"

ENERGY_TEXT="Ten Tau (200 PU) | 13.6 TeV"

DELTAR_DIRS=(
    "DeltaR0p3"
    "DeltaR0p25"
    "DeltaR0p2"
    "DeltaR0p15"
    "DeltaR0p1"
)

DELTAR_LABELS=(
    '$\Delta R=0.3$'
    '$\Delta R=0.25$'
    '$\Delta R=0.2$'
    '$\Delta R=0.15$'
    '$\Delta R=0.1$'
)

PT_REBIN="0,5,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,220,240,260,280,300"
ETA_REBIN="-2.4,-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0,2.4"
PHI_REBIN="-3.5,-2.8,-2.1,-1.4,-0.7,0.0,0.7,1.4,2.1,2.8,3.5"

has_hist() {
    local hist="$1"
    rootls "$DQM_FILE:$hist" >/dev/null 2>&1
}

join_by_comma() {
    local IFS=","
    echo "$*"
}

run_cmd() {
    echo
    echo "Running:"
    echo "$*"
    "$@"
}

make_comparison_plot() {
    local hist_base="$1"
    local name="$2"
    local xlabel="$3"
    local ylabel="$4"
    local xlim="$5"
    local ylim="$6"
    local ylim_ratio="$7"
    local rebin="$8"
    local inverted="${9:-0}"

    local files=()
    local hists=()
    local labels=()

    for i in "${!DELTAR_DIRS[@]}"; do
        local hist="${BASE_DIR}/${DELTAR_DIRS[$i]}/${hist_base}"

        if has_hist "$hist"; then
            files+=("$DQM_FILE")
            hists+=("$hist")
            labels+=("${DELTAR_LABELS[$i]}")
        else
            echo "Skipping missing histogram: $hist"
        fi
    done

    if [ "${#hists[@]}" -eq 0 ]; then
        echo "No valid histograms for ${hist_base}. Skipping."
        return
    fi

    local cmd=(
        python3 "$MAKE_COMPARISON"
        --files "$(join_by_comma "${files[@]}")"
        --hists "$(join_by_comma "${hists[@]}")"
        --labels "$(join_by_comma "${labels[@]}")"
        --xlabel "$xlabel"
        --ylabel "$ylabel"
        --energy-text "$ENERGY_TEXT"
        --odir "$OUTDIR_COMPARISON"
        --name "$name"
    )

    if [ -n "$xlim" ]; then
        cmd+=("--xlim=${xlim}")
    fi

    if [ -n "$ylim" ]; then
        cmd+=("--ylim=${ylim}")
    fi

    if [ -n "$ylim_ratio" ]; then
        cmd+=("--ylim-ratio=${ylim_ratio}")
    fi

    if [ -n "$rebin" ]; then
        cmd+=("--rebin=${rebin}")
    fi

    if [ "$inverted" -eq 1 ]; then
        cmd+=(--inverted)
    fi

    run_cmd "${cmd[@]}"
}

make_summary_plot() {
    local delta_dir="$1"
    local delta_text="$2"
    local variable="$3"
    local den_hist="$4"
    local num_hist="$5"
    local rate_hist="$6"
    local name="$7"
    local xlabel="$8"
    local ylabel="$9"
    local xlim="${10}"
    local ylim="${11}"
    local rebin="${12}"
    local rate_label="${13}"
    local den_label="${14}"
    local num_label="${15}"
    local extra_args="${16:-}"

    local den="${BASE_DIR}/${delta_dir}/${den_hist}"
    local num="${BASE_DIR}/${delta_dir}/${num_hist}"
    local rate="${BASE_DIR}/${delta_dir}/${rate_hist}"

    if ! has_hist "$den"; then
        echo "Skipping missing denominator: $den"
        return
    fi

    if ! has_hist "$num"; then
        echo "Skipping missing numerator: $num"
        return
    fi

    if ! has_hist "$rate"; then
        echo "Skipping missing rate: $rate"
        return
    fi

    local cmd=(
        python3 "$MAKE_TAU_VALIDATION"
        --mode summary
        --files "$DQM_FILE"
        --den-hists "$den"
        --num-hists "$num"
        --rate-hists "$rate"
        --labels "$rate_label"
        --den-label "$den_label"
        --num-label "$num_label"
        --xlabel "$xlabel"
        --ylabel "$ylabel"
        --ylim "$ylim"
        --energy-text "$ENERGY_TEXT"
        --text "Tau validation ${delta_text}"
        --odir "${OUTDIR_SUMMARY}/${delta_dir}"
        --name "$name"
    )

    if [ -n "$xlim" ]; then
        cmd+=("--xlim=${xlim}")
    fi

    if [ -n "$rebin" ]; then
        cmd+=("--rebin=${rebin}")
    fi

    if [ -n "$extra_args" ]; then
        cmd+=("$extra_args")
    fi

    run_cmd "${cmd[@]}"
}

make_response_plot() {
    local base="$1"
    local name="$2"
    local xlabel="$3"
    local ylabel="$4"
    local xlim="$5"
    local ylim="$6"
    local ylim_ratio="$7"
    local rebin="$8"

    local files=()
    local mean_hists=()
    local sigma_hists=()
    local labels=()

    for i in "${!DELTAR_DIRS[@]}"; do
        local mean_hist="${BASE_DIR}/${DELTAR_DIRS[$i]}/${base}_Mean"
        local sigma_hist="${BASE_DIR}/${DELTAR_DIRS[$i]}/${base}_Sigma"

        if has_hist "$mean_hist" && has_hist "$sigma_hist"; then
            files+=("$DQM_FILE")
            mean_hists+=("$mean_hist")
            sigma_hists+=("$sigma_hist")
            labels+=("${DELTAR_LABELS[$i]}")
        else
            echo "Skipping missing response pair:"
            echo "  mean : $mean_hist"
            echo "  sigma: $sigma_hist"
        fi
    done

    if [ "${#mean_hists[@]}" -eq 0 ]; then
        echo "No valid response histograms for ${base}. Skipping."
        return
    fi

    local cmd=(
        python3 "$MAKE_TAU_VALIDATION"
        --mode response
        --files "$(join_by_comma "${files[@]}")"
        --mean-hists "$(join_by_comma "${mean_hists[@]}")"
        --sigma-hists "$(join_by_comma "${sigma_hists[@]}")"
        --labels "$(join_by_comma "${labels[@]}")"
        --xlabel "$xlabel"
        --ylabel "$ylabel"
        --energy-text "$ENERGY_TEXT"
        --odir "$OUTDIR_RESPONSE"
        --name "$name"
        --title "$name"
    )

    if [ -n "$xlim" ]; then
        cmd+=("--xlim=${xlim}")
    fi

    if [ -n "$ylim" ]; then
        cmd+=("--ylim=${ylim}")
    fi

    if [ -n "$ylim_ratio" ]; then
        cmd+=("--ylim-ratio=${ylim_ratio}")
    fi

    if [ -n "$rebin" ]; then
        cmd+=("--rebin=${rebin}")
    fi

    run_cmd "${cmd[@]}"
}

mkdir -p "$OUTDIR_COMPARISON"
mkdir -p "$OUTDIR_SUMMARY"
mkdir -p "$OUTDIR_RESPONSE"

echo "Input file:"
echo "$DQM_FILE"

echo
echo "Making comparison plots"

make_comparison_plot "Eff_vs_pt" "Eff_vs_pt_DeltaR_comparison" '$p_T$ [GeV]' "Efficiency" "0,300" "" "" "$PT_REBIN"
make_comparison_plot "Eff_vs_eta" "Eff_vs_eta_DeltaR_comparison" '$\eta$' "Efficiency" "-2.5,2.5" "" "" "$ETA_REBIN"
make_comparison_plot "Eff_vs_phi" "Eff_vs_phi_DeltaR_comparison" '$\phi$' "Efficiency" "" "" "" "$PHI_REBIN"

make_comparison_plot "Fake_vs_pt" "Fake_vs_pt_DeltaR_comparison" '$p_T$ [GeV]' "Fake rate" "0,300" "" "" "$PT_REBIN" 1
make_comparison_plot "Fake_vs_eta" "Fake_vs_eta_DeltaR_comparison" '$\eta$' "Fake rate" "-2.5,2.5" "" "" "$ETA_REBIN" 1
make_comparison_plot "Fake_vs_phi" "Fake_vs_phi_DeltaR_comparison" '$\phi$' "Fake rate" "" "" "" "$PHI_REBIN" 1

make_comparison_plot "Dup_vs_pt" "Dup_vs_pt_DeltaR_comparison" '$p_T$ [GeV]' "Duplicate rate" "0,300" "" "" "$PT_REBIN"
make_comparison_plot "Dup_vs_eta" "Dup_vs_eta_DeltaR_comparison" '$\eta$' "Duplicate rate" "-2.5,2.5" "" "" "$ETA_REBIN"
make_comparison_plot "Dup_vs_phi" "Dup_vs_phi_DeltaR_comparison" '$\phi$' "Duplicate rate" "" "" "" "$PHI_REBIN"

make_comparison_plot "Split_vs_pt" "Split_vs_pt_DeltaR_comparison" '$p_T$ [GeV]' "Split rate" "0,300" "" "" "$PT_REBIN"
make_comparison_plot "Split_vs_eta" "Split_vs_eta_DeltaR_comparison" '$\eta$' "Split rate" "-2.5,2.5" "" "" "$ETA_REBIN"
make_comparison_plot "Split_vs_phi" "Split_vs_phi_DeltaR_comparison" '$\phi$' "Split rate" "" "" "" "$PHI_REBIN"

echo
echo "Making summary plots"

for i in "${!DELTAR_DIRS[@]}"; do
    delta_dir="${DELTAR_DIRS[$i]}"
    delta_text="${DELTAR_LABELS[$i]}"

    make_summary_plot "$delta_dir" "$delta_text" "pt" "genTau_pt" "genTauMatched_pt" "Eff_vs_pt" "Tau_Efficiency_pt" '$p_T$ [GeV]' "Tau Efficiency" "0,300" "0,1.3" "$PT_REBIN" "Efficiency" "Gen taus" "Gen taus matched to reco taus"
    make_summary_plot "$delta_dir" "$delta_text" "eta" "genTau_eta" "genTauMatched_eta" "Eff_vs_eta" "Tau_Efficiency_eta" '$\eta$' "Tau Efficiency" "-2.5,2.5" "0,1.3" "$ETA_REBIN" "Efficiency" "Gen taus" "Gen taus matched to reco taus"
    make_summary_plot "$delta_dir" "$delta_text" "phi" "genTau_phi" "genTauMatched_phi" "Eff_vs_phi" "Tau_Efficiency_phi" '$\phi$' "Tau Efficiency" "" "0,1.3" "$PHI_REBIN" "Efficiency" "Gen taus" "Gen taus matched to reco taus"

    make_summary_plot "$delta_dir" "$delta_text" "pt" "recoTau_pt" "recoTauMatched_pt" "Fake_vs_pt" "Tau_FakeRate_pt" '$p_T$ [GeV]' "Tau Fake Rate" "0,300" "0,1.3" "$PT_REBIN" "Fake rate" "Reco taus" "Reco taus matched to gen taus" "--inverted"
    make_summary_plot "$delta_dir" "$delta_text" "eta" "recoTau_eta" "recoTauMatched_eta" "Fake_vs_eta" "Tau_FakeRate_eta" '$\eta$' "Tau Fake Rate" "-2.5,2.5" "0,1.3" "$ETA_REBIN" "Fake rate" "Reco taus" "Reco taus matched to gen taus" "--inverted"
    make_summary_plot "$delta_dir" "$delta_text" "phi" "recoTau_phi" "recoTauMatched_phi" "Fake_vs_phi" "Tau_FakeRate_phi" '$\phi$' "Tau Fake Rate" "" "0,1.3" "$PHI_REBIN" "Fake rate" "Reco taus" "Reco taus matched to gen taus" "--inverted"
done

echo
echo "Making response plots"

make_response_plot "ResponsePt_RecoOverGen_vs_pt" "ResponsePt_resolution_vs_pt_DeltaR_comparison" '$p_T^{gen}$ [GeV]' '$\sigma(p_T^{reco}/p_T^{gen}) / \langle p_T^{reco}/p_T^{gen} \rangle$' "0,300" "" "" "$PT_REBIN"
make_response_plot "ResponsePt_RecoOverGen_vs_eta" "ResponsePt_resolution_vs_eta_DeltaR_comparison" '$\eta^{gen}$' '$\sigma(p_T^{reco}/p_T^{gen}) / \langle p_T^{reco}/p_T^{gen} \rangle$' "-2.5,2.5" "" "" "$ETA_REBIN"
make_response_plot "ResponsePt_RecoOverGen_vs_phi" "ResponsePt_resolution_vs_phi_DeltaR_comparison" '$\phi^{gen}$' '$\sigma(p_T^{reco}/p_T^{gen}) / \langle p_T^{reco}/p_T^{gen} \rangle$' "" "" "" ""
make_response_plot "ResponsePt_RecoOverGen_vs_mass" "ResponsePt_resolution_vs_mass_DeltaR_comparison" '$m^{gen}$ [GeV]' '$\sigma(p_T^{reco}/p_T^{gen}) / \langle p_T^{reco}/p_T^{gen} \rangle$' "" "" "" "$PHI_REBIN"

make_response_plot "ResponseMass_RecoOverGen_vs_pt" "ResponseMass_resolution_vs_pt_DeltaR_comparison" '$p_T^{gen}$ [GeV]' '$\sigma(m^{reco}/m^{gen}) / \langle m^{reco}/m^{gen} \rangle$' "0,300" "" "" "$PT_REBIN"
make_response_plot "ResponseMass_RecoOverGen_vs_eta" "ResponseMass_resolution_vs_eta_DeltaR_comparison" '$\eta^{gen}$' '$\sigma(m^{reco}/m^{gen}) / \langle m^{reco}/m^{gen} \rangle$' "-2.5,2.5" "" "" "$ETA_REBIN"
make_response_plot "ResponseMass_RecoOverGen_vs_phi" "ResponseMass_resolution_vs_phi_DeltaR_comparison" '$\phi^{gen}$' '$\sigma(m^{reco}/m^{gen}) / \langle m^{reco}/m^{gen} \rangle$' "" "" "" "$PHI_REBIN"
make_response_plot "ResponseMass_RecoOverGen_vs_mass" "ResponseMass_resolution_vs_mass_DeltaR_comparison" '$m^{gen}$ [GeV]' '$\sigma(m^{reco}/m^{gen}) / \langle m^{reco}/m^{gen} \rangle$' "" "" "" "$PHI_REBIN"

echo
echo "Done."