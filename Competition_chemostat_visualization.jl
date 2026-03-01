# Reads summary CSV with classification (string) and consensus
# Supports multi-subplot layout for different D values or parameter variations

using CSV, DataFrames, Plots, Printf
using Plots.PlotMeasures
using ColorSchemes
using LaTeXStrings

# ==============================================================================
# Classification String → Numeric Code Mapping
# ==============================================================================

"""
    classification_to_code(classification::String) -> Int

Convert classification string from summary CSV back to numeric winner code
for color mapping compatibility.

Mapping:
  "0 win"                → 0
  "1 win"                → 1
  "2 win"                → 2
  "901 coexist"          → 901   (strategy 0 & 1 coexist, 0 ≥ 1)
  "10 coexist"           → 10    (strategy 0 & 1 coexist, 1 > 0)
  "12 coexist"           → 12    (strategy 1 & 2 coexist, 1 > 2)
  "21 coexist"           → 21    (strategy 1 & 2 coexist, 2 > 1)
  "02 coexist"           → 02    (strategy 0 & 2 coexist, 0 > 2)
  "20 coexist"           → 20    (strategy 0 & 2 coexist, 2 > 0)
  "Extinction"           → -1
  "Neutral coexistence"  → 99
  "Multi-stability"      → 98
  "No data" / unknown    → -999
"""
function classification_to_code(classification::AbstractString)
    cl = strip(classification)

    # Exclusion: "<number> win"
    m = match(r"^(\d+)\s+win$", cl)
    if m !== nothing
        return parse(Int, m.captures[1])
    end

    # Coexistence: "<number> coexist"
    m = match(r"^(\d+)\s+coexist$", cl)
    if m !== nothing
        return parse(Int, m.captures[1])
    end

    cl == "Extinction"          && return -1
    cl == "Neutral coexistence" && return 99
    cl == "Multi-stability"     && return 98
    cl == "No data"             && return -999

    return -999  # fallback
end

# ==============================================================================
# Data Loading
# ==============================================================================

"""
Load data from multiple summary CSV files.
Each file becomes one subplot.
"""
function load_multiple_csv_files(csv_files)
    println("=== Loading Summary Data from Multiple Files ===")

    file_data_list = []

    for (i, csv_file) in enumerate(csv_files)
        if !isfile(csv_file)
            println("Warning: CSV file not found: $(csv_file)")
            continue
        end

        println("Loading file $(i)/$(length(csv_files)): $(basename(csv_file))")

        df = CSV.read(csv_file, DataFrame;
                      stringtype=String,
                      dateformat=nothing,
                      types=Dict("strategy_pair" => String))

        # Get D values
        if "D" in names(df)
            file_D_values = sort(unique(df.D))
        elseif "filename_D" in names(df)
            file_D_values = sort(unique(df.filename_D))
        else
            println("Skipping file: no D column found")
            continue
        end

        println("  - Loaded $(nrow(df)) records with D values: $(file_D_values)")

        push!(file_data_list, Dict(
            "df" => df, "file_index" => i, "file_name" => basename(csv_file),
            "D_values" => file_D_values
        ))
    end

    if isempty(file_data_list); error("No valid CSV files found!"); end
    println("Total files processed: $(length(file_data_list))")
    return organize_file_data(file_data_list)
end

"""
Organize file data into grids for plotting.
Handles both old format (dominant_winner/confidence) and new format (classification/consensus).
"""
function organize_file_data(file_data_list)
    # Clean strategy_pair column (fix date-like parsing artifacts)
    for file_data in file_data_list
        df = file_data["df"]
        if "strategy_pair" in names(df)
            df.strategy_pair = string.(df.strategy_pair)
            df.strategy_pair = replace.(df.strategy_pair, r"\d{4}-01-01" => "0-1")
            df.strategy_pair = replace.(df.strategy_pair, r"\d{4}-02-01" => "1-2")
            df.strategy_pair = replace.(df.strategy_pair, r"\d{4}-03-01" => "0-2")
        end
    end

    # Collect all T1, T2 values across files
    all_T1_values = Set(); all_T2_values = Set()
    for file_data in file_data_list
        df = file_data["df"]
        if "T1" in names(df) && "T2" in names(df)
            union!(all_T1_values, df.T1); union!(all_T2_values, df.T2)
        elseif "filename_T1" in names(df) && "filename_T2" in names(df)
            union!(all_T1_values, df.filename_T1); union!(all_T2_values, df.filename_T2)
        end
    end

    T1_values = sort(collect(all_T1_values))
    T2_values = sort(collect(all_T2_values))
    println("Combined T1 range: $(T1_values[1]) to $(T1_values[end]) ($(length(T1_values)) values)")
    println("Combined T2 range: $(T2_values[1]) to $(T2_values[end]) ($(length(T2_values)) values)")

    winner_grids = Dict()
    file_identifiers = []

    for file_data in file_data_list
        df = file_data["df"]
        file_idx = file_data["file_index"]
        file_id = "file_$(file_idx)"
        push!(file_identifiers, file_id)

        # Determine column names
        T1_col, T2_col = if "T1" in names(df); (:T1, :T2) else (:filename_T1, :filename_T2) end

        # Detect format: new (classification string) or old (numeric winner)
        use_classification = "classification" in names(df)

        winner_col = if use_classification
            :classification
        elseif "dominant_winner" in names(df)
            :dominant_winner
        elseif "dominant_winner_code" in names(df)
            :dominant_winner_code
        else
            :winner
        end

        winner_grid = fill(-999, length(T1_values), length(T2_values))

        for row in eachrow(df)
            t1_idx = findfirst(x -> x == row[T1_col], T1_values)
            t2_idx = findfirst(x -> x == row[T2_col], T2_values)
            if t1_idx !== nothing && t2_idx !== nothing
                if use_classification
                    # New format: string → numeric code
                    winner_grid[t1_idx, t2_idx] = classification_to_code(string(row[winner_col]))
                else
                    # Old format: numeric code directly
                    original_winner = row[winner_col]

                    # T1=0 coexistence → Neutral coexist
                    processed_winner = if row[T1_col] == 0 && is_coexistence_code(original_winner)
                        99
                    else
                        original_winner
                    end

                    winner_grid[t1_idx, t2_idx] = processed_winner
                end
            end
        end

        winner_grids[file_id] = winner_grid
        println("File $(file_idx): filled $(count(x -> x != -999, winner_grid)) / $(length(winner_grid)) grid points")
    end

    # Detect strategy_pair
    strategy_pair = "unknown"
    if !isempty(file_data_list)
        first_df = file_data_list[1]["df"]
        if "strategy_pair" in names(first_df)
            unique_pairs = unique(first_df.strategy_pair)
            if length(unique_pairs) == 1; strategy_pair = unique_pairs[1]; end
            println("Found strategy_pair: $(strategy_pair)")
        end
    end

    return Dict(
        "T1_values" => T1_values, "T2_values" => T2_values,
        "D_values" => file_identifiers,
        "winner_grids" => winner_grids,
        "strategy_pair" => strategy_pair,
        "file_data_list" => file_data_list
    )
end

"""
Check if a numeric winner code represents coexistence (for old format compatibility).
"""
function is_coexistence_code(winner)
    coexistence_codes = [10, 12, 21, 20, 901]
    if winner in coexistence_codes
        return true
    end
    # Three-digit codes (e.g. 102) are coexistence
    if winner >= 100 && winner <= 999
        return true
    end
    return false
end

# ==============================================================================
# Color Mapping & Labels
# ==============================================================================

"""
Interpret winner code to display label.
Strategy 0 = N, Strategy 1 = T, Strategy 2 = P
"""
function interpret_winner_code(winner, strategy_pair="unknown")
    winner == 99   && return "Neutral coexistence"
    winner == 98   && return "Multi-stability"
    winner == -1   && return "Extinct"
    winner == -999 && return "No data"
    winner == 0    && return "N wins"
    winner == 1    && return "T wins"
    winner == 2    && return "P wins"
    winner == 901  && return "Stable coex. (N > T)"
    winner == 10   && return "Stable coex. (T > N)"
    winner == 12   && return "Stable coex. (T > P)"
    winner == 21   && return "Stable coex. (P > T)"
    winner == 20   && return "Stable coex. (P > N)"
    return "Other ($(winner))"
end

"""
Create color mapping for winner codes.
Returns Dict(winner_code => color_index).
"""
function get_color_mapping(strategy_pair="unknown")
    color_map = Dict(
        -999 => 1,   # Missing/Failed → light grey
        -998 => 1,
        -997 => 1,
        -1   => 2,   # Extinction → dark grey
        99   => 11,  # Neutral coexist → grey-blue
        98   => 12,  # Multi-stability → amber/brown
    )

    if strategy_pair == "0-1" || strategy_pair == "01"
        color_map[0]   = 7   # N → orange
        color_map[1]   = 8   # T → cyan
        color_map[901] = 4   # N&T coexist (N≥T) → purple-pink
        color_map[10]  = 3   # T&N coexist (T>N) → blue-purple
    elseif strategy_pair == "1-2" || strategy_pair == "12"
        color_map[1]  = 8    # T → cyan
        color_map[2]  = 9    # P → salmon
        color_map[12] = 3    # T&P coexist (T>P) → blue-purple
        color_map[21] = 4    # P&T coexist (P>T) → purple-pink
    elseif strategy_pair == "0-2" || strategy_pair == "02"
        color_map[0]  = 7    # N → orange
        color_map[2]  = 9    # P → salmon
        color_map[20] = 3    # P&N coexist (P>N) → blue-purple
    else
        # Generic fallback
        color_map[0] = 7; color_map[1] = 8; color_map[2] = 9
        color_map[901] = 4; color_map[10] = 3
        color_map[12] = 3; color_map[21] = 4
        color_map[20] = 3
    end

    return color_map
end

"""
Color palette (12 colors, indexed 1-12).
"""
function get_color_palette()
    return [
        RGB(0.9, 0.9, 0.9),    #  1: Missing/Failed (light grey)
        RGB(0.5, 0.5, 0.5),    #  2: Extinction (dark grey)
        RGB(0.4, 0.6, 0.8),    #  3: Coexist variant A (blue-purple)
        RGB(0.8, 0.6, 0.8),    #  4: Coexist variant B (purple-pink)
        RGB(0.8, 0.6, 0.8),    #  5: (unused)
        RGB(0.85, 0.8, 0.95),  #  6: (unused)
        RGB(0.95, 0.8, 0.25),  #  7: Strategy N wins (orange-yellow)
        RGB(0.35, 0.8, 0.9),   #  8: Strategy T wins (cyan)
        RGB(0.95, 0.6, 0.5),   #  9: Strategy P wins (salmon)
        RGB(0.8, 0.4, 0.6),    # 10: (unused)
        RGB(0.6, 0.7, 0.75),   # 11: Neutral coexist (grey-blue)
        RGB(0.75, 0.55, 0.2),  # 12: Multi-stability (amber/brown)
    ]
end

# ==============================================================================
# Helper Functions
# ==============================================================================

function auto_detect_intervals(parameter_values)
    if length(parameter_values) < 2
        return 1.0
    end
    interval = parameter_values[2] - parameter_values[1]
    return interval
end

function smart_tick_selection(values; max_ticks=8, min_ticks=3)
    n_values = length(values)

    if n_values <= max_ticks
        return values, collect(1:n_values)
    end

    step = max(1, div(n_values, max_ticks))
    selected_indices = [1]

    current_idx = 1 + step
    while current_idx < n_values
        push!(selected_indices, current_idx)
        current_idx += step
    end

    if selected_indices[end] != n_values
        push!(selected_indices, n_values)
    end

    if length(selected_indices) < min_ticks && n_values >= min_ticks
        step = max(1, div(n_values, min_ticks))
        selected_indices = collect(1:step:n_values)
        if selected_indices[end] != n_values
            push!(selected_indices, n_values)
        end
    end

    selected_values = values[selected_indices]
    return selected_values, selected_indices
end

# ==============================================================================
# Subplot Label Letters
# ==============================================================================

const SUBPLOT_LABELS = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)"]

"""
Add subplot label (e.g. "(a)") outside the plot area, at the top-left corner.
"""
# for 8
#function add_subplot_label!(p, label; fontsize=22, x_offset=-0.1, y_offset=1.05)
#for 4
#function add_subplot_label!(p, label; fontsize=22, x_offset=-0.35, y_offset=1.08)
#for 1
#function add_subplot_label!(p, label; fontsize=22, x_offset=-0.4, y_offset=1.08)
#for 6
function add_subplot_label!(p, label; fontsize=22, x_offset=-0.35, y_offset=1.08)
    xl = xlims(p)
    yl = ylims(p)
    x_pos = xl[1] + x_offset * (xl[2] - xl[1])
    y_pos = yl[1] + y_offset * (yl[2] - yl[1])
    annotate!(p, [(x_pos, y_pos, text(label, :left, :bottom, fontsize, :bold))])
end

# ==============================================================================
# Legend Creation
# ==============================================================================

"""
Define the fixed display order of legend items for each strategy pair.
Row 1: exclusive winners + coexistence states;  Row 2: extinction + neutral + multi-stability.
Returns Vector{Vector{Int}} — each inner vector is one row of winner codes.
"""
function get_legend_display_order(strategy_pair)
    if strategy_pair == "0-1" || strategy_pair == "01"
        return [[1, 0, 901, 10],    # T wins, N wins, Stable coex. (N>T), Stable coex. (T>N)
                [-1, 99, 98]]       # Extinct, Neutral coexistence, Multi-stability
    elseif strategy_pair == "1-2" || strategy_pair == "12"
        return [[1, 2, 12, 21],     # T wins, P wins, Stable coex. (T>P), Stable coex. (P>T)
                [-1, 99, 98]]       # Extinct, Neutral coexistence, Multi-stability
    elseif strategy_pair == "0-2" || strategy_pair == "02"
        return [[0, 2, 20],         # N wins, P wins, Stable coex. (P>N)
                [-1, 99, 98]]       # Extinct, Neutral coexistence, Multi-stability
    else
        # Generic fallback
        return [[0, 1, 2, 901, 10, 12, 21, 20],
                [-1, 99, 98]]
    end
end

function create_legend(strategy_pair, color_map, colors, cmap, winner_grids, D_values;
                      legend_fontsize=18, rect_scale=0.6, legend_width=1772, legend_height=300,
                      max_per_row=4)
    # Collect all winner codes actually present in the data
    all_winners = Set()
    for D in D_values
        winner_grid = winner_grids[D]
        unique_winners = unique(vec(winner_grid))
        union!(all_winners, filter(x -> x != -999, unique_winners))
    end

    # Get fixed display order, then filter to only codes present in data
    display_rows = get_legend_display_order(strategy_pair)
    filtered_rows = [filter(w -> w in all_winners && haskey(color_map, w), row) for row in display_rows]
    filtered_rows = filter(!isempty, filtered_rows)   # remove empty rows

    # Smart layout: if total items ≤ max_per_row, collapse all into a single row
    total_items = sum(length.(filtered_rows))
    if total_items <= max_per_row
        filtered_rows = [reduce(vcat, filtered_rows)]
    end

    n_rows = length(filtered_rows)
    if n_rows == 0
        return plot(xlims=(0,1), ylims=(0,1), framestyle=:none, grid=false, legend=false,
                    size=(legend_width, legend_height))
    end

    legend_plot = plot(xlims=(0,1), ylims=(0,1),
                      framestyle=:none, grid=false, legend=false,
                      size=(legend_width, legend_height),
                      margin=5Plots.mm)

    # Column positions: single row → evenly distribute by actual count;
    #                   multi-row  → fixed grid aligned to max_per_row
    rect_height = n_rows == 1 ? 0.20 : 0.14
    n_cols_for_pos = n_rows == 1 ? length(filtered_rows[1]) : max_per_row
    col_x_positions = if n_rows == 1 && n_cols_for_pos > 1
        range(0.20, 0.80, length=n_cols_for_pos)   # narrower range to avoid label clipping
    elseif n_cols_for_pos == 1
        [0.50]                                       # single item centered
    else
        range(0.10, 0.90, length=n_cols_for_pos)     # multi-row: full width for alignment
    end
    rect_width  = 0.8 / n_cols_for_pos * rect_scale

    for (row_idx, row_codes) in enumerate(filtered_rows)
        # Y positions for this row
        if n_rows == 1
            y_box   = 0.60
            y_label = 0.22
        else
            y_box   = 1.0 - row_idx * (0.85 / n_rows) + 0.10
            y_label = y_box - rect_height / 2 - 0.10
        end

        for (col_idx, winner_code) in enumerate(row_codes)
            if col_idx > n_cols_for_pos; break; end

            x_center = col_x_positions[col_idx]
            color_idx = color_map[winner_code]
            colormap_value = (color_idx - 0.5) / length(colors)
            exact_color = cmap[colormap_value]

            # Draw color box
            plot!(legend_plot,
                  Shape([x_center - rect_width/2, x_center + rect_width/2,
                         x_center + rect_width/2, x_center - rect_width/2],
                        [y_box - rect_height/2, y_box - rect_height/2,
                         y_box + rect_height/2, y_box + rect_height/2]),
                  color=exact_color, linecolor=:black, linewidth=1)

            # Draw label
            label_text = interpret_winner_code(winner_code, strategy_pair)
            annotate!(legend_plot, x_center, y_label,
                     text(label_text, :center, legend_fontsize))
        end
    end

    return legend_plot
end

# ==============================================================================
# Main Plotting Function
# ==============================================================================

function plot_robustness_multi_subplots(results;
                                      save_path="robustness_phase_diagram.pdf",
                                      legend_fontsize=18,
                                      rect_scale=0.6,
                                      max_x_ticks=4,
                                      max_y_ticks=4,
                                      subplot_title_fontsize=22,
                                      axis_label_fontsize=20,
                                      tick_fontsize=18,
                                      # 150mm @ 300dpi ≈ 1772 pixels (matching LaTeX textwidth)
                                      figure_width=1772,
                                      legend_height=300,
                                      subplot_label_fontsize=22,
                                      show_subplot_labels=true,
                                      title_gap=0,
                                      dpi=300,
                                      max_plots=9,
                                      # Custom subplot titles: pass a Vector{String} to override auto titles
                                      # e.g. subplot_titles=["Low D", "Medium D", "High D"]
                                      # If nothing (default), auto-generates "D = x.xxx" from data
                                      subplot_titles=nothing)

    T1_values = results["T1_values"]
    T2_values = results["T2_values"]
    D_values = results["D_values"]
    strategy_pair = results["strategy_pair"]

    D_subset = D_values[1:min(max_plots, length(D_values))]
    n_plots = length(D_subset)

    # Determine layout
    if n_plots == 1
        layout_tuple = (1, 1)
    elseif n_plots == 2
        layout_tuple = (1, 2)
    elseif n_plots <= 3
        layout_tuple = (1, 3)
    elseif n_plots <= 4
        layout_tuple = (1, 4)
    elseif n_plots <= 6
        layout_tuple = (2, 3)
    elseif n_plots <= 8
        layout_tuple = (2, 4)
    else
        layout_tuple = (3, 3)
    end

    # Calculate figure dimensions and font sizes
    n_rows, n_cols = layout_tuple

    if n_plots == 1
        # Single plot: fixed size and font sizes tuned for small figure
        effective_width = 500
        total_height    = 600
        subplot_frac    = 0.80
        s_title_fs  = 12;  s_axis_fs  = 10;  s_tick_fs    = 9
        s_legend_fs = 9;   s_label_fs = 12
    else
        # Multi-plot: reference-based sizing (2-row case as baseline)
        effective_width  = figure_width
        ref_row_height   = round(Int, figure_width / n_cols * 0.85)
        ref_total_height = 2 * ref_row_height + legend_height
        ref_per_row_px   = 0.85 * ref_total_height / 2
        ref_legend_px    = 0.15 * ref_total_height
        subplot_area_px  = n_rows * ref_per_row_px
        total_height     = round(Int, subplot_area_px + ref_legend_px)
        subplot_frac     = subplot_area_px / (subplot_area_px + ref_legend_px)
        s_title_fs  = subplot_title_fontsize;  s_axis_fs  = axis_label_fontsize
        s_tick_fs   = tick_fontsize;            s_legend_fs = legend_fontsize
        s_label_fs  = subplot_label_fontsize
    end

    # Generate title gap string
    gap_str = repeat("\n", title_gap)

    println("=== Creating Phase Diagram Multi-Subplot Plot ===")
    println("Strategy pair: $(strategy_pair)")
    println("Files: $(length(D_subset))")
    println("Layout: $(layout_tuple[1])×$(layout_tuple[2]) for $(n_plots) subplots")

    # Color setup
    color_map = get_color_mapping(strategy_pair)
    colors = get_color_palette()
    cmap = cgrad(colors; categorical=true)

    subplots = []

    for (idx, D) in enumerate(D_subset)
        println("Creating subplot $(idx)/$(length(D_subset)) for $(D)")

        # Get file-specific data range
        file_info = results["file_data_list"][findfirst(
            fd -> fd["file_index"] == parse(Int, split(string(D), "_")[2]),
            results["file_data_list"])]
        file_df = file_info["df"]

        T1_col = "T1" in names(file_df) ? :T1 : :filename_T1
        T2_col = "T2" in names(file_df) ? :T2 : :filename_T2

        file_T1_values = sort(unique(file_df[!, T1_col]))
        file_T2_values = sort(unique(file_df[!, T2_col]))

        # Boundaries for heatmap cells
        file_T1_interval = auto_detect_intervals(file_T1_values)
        file_T2_interval = auto_detect_intervals(file_T2_values)

        file_T1_boundaries = [file_T1_values[1] - file_T1_interval/2;
                              [v + file_T1_interval/2 for v in file_T1_values]]
        file_T2_boundaries = [file_T2_values[1] - file_T2_interval/2;
                              [v + file_T2_interval/2 for v in file_T2_values]]

        # Tick selection
        selected_file_T1_values, _ = smart_tick_selection(file_T1_values; max_ticks=max_y_ticks)
        selected_file_T2_values, _ = smart_tick_selection(file_T2_values; max_ticks=max_x_ticks)

        # Build heatmap data
        grid_data = results["winner_grids"][D]
        t1_indices = [findfirst(x -> x == t1, T1_values) for t1 in file_T1_values]
        t2_indices = [findfirst(x -> x == t2, T2_values) for t2 in file_T2_values]

        heat_data = zeros(length(file_T1_values), length(file_T2_values))
        for i in 1:length(file_T1_values), j in 1:length(file_T2_values)
            winner = grid_data[t1_indices[i], t2_indices[j]]
            heat_data[i, j] = get(color_map, winner, 1)
        end

        # Subplot title: use custom title if provided, otherwise auto-generate
        title_text = if subplot_titles !== nothing && idx <= length(subplot_titles)
            subplot_titles[idx]
        else
            @sprintf("D = %.3f", file_info["D_values"][1])
        end

        # Determine if this subplot should show axis labels
        show_xlabel = (idx-1) >= (layout_tuple[1]-1) * layout_tuple[2]
        show_ylabel = (idx-1) % layout_tuple[2] == 0

        # Create base plot with axes and styling (no heatmap — use Shape for crisp PDF)
        p = plot(
            title=title_text * gap_str,
            xlabel=show_xlabel ? "T₂ (S₁,ᵢₙ OFF Duration)" : "",
            ylabel=show_ylabel ? "T₁ (S₁,ᵢₙ ON Duration)" : "",
            framestyle=:box,
            xticks=(selected_file_T2_values,
                    [isinteger(v) ? string(Int(v)) : string(v) for v in selected_file_T2_values]),
            yticks=selected_file_T1_values,
            xlims=(file_T2_boundaries[1], file_T2_boundaries[end]),
            ylims=(file_T1_boundaries[1], file_T1_boundaries[end]),
            legend=false,
            titlefontsize=s_title_fs,
            xlabelfontsize=s_axis_fs, ylabelfontsize=s_axis_fs,
            xtickfontsize=s_tick_fs, ytickfontsize=s_tick_fs
        )

        # Draw each cell as a solid-colored rectangle (vector graphics, no interpolation)
        n_T1, n_T2 = size(heat_data)
        for j in 1:n_T2
            for i in 1:n_T1
                c = colors[clamp(Int(heat_data[i, j]), 1, length(colors))]
                plot!(p, Shape(
                    [file_T2_boundaries[j], file_T2_boundaries[j+1],
                     file_T2_boundaries[j+1], file_T2_boundaries[j]],
                    [file_T1_boundaries[i], file_T1_boundaries[i],
                     file_T1_boundaries[i+1], file_T1_boundaries[i+1]]
                ), fillcolor=c, linecolor=c, linewidth=0, label=false)
            end
        end

        # Add (a), (b), (c) ... label outside the plot (top-left)
        if show_subplot_labels && idx <= length(SUBPLOT_LABELS)
            add_subplot_label!(p, SUBPLOT_LABELS[idx]; fontsize=s_label_fs)
        end

        push!(subplots, p)
    end

    # Fill remaining layout positions with empty plots
    total_positions = layout_tuple[1] * layout_tuple[2]
    while length(subplots) < total_positions
        push!(subplots, plot(framestyle=:none, grid=false, legend=false))
    end

    # Assemble subplots (no main title)
    subplot_layout = plot(subplots...,
                         layout=layout_tuple,
                         left_margin=12Plots.mm, right_margin=8Plots.mm,
                         bottom_margin=10Plots.mm, top_margin=8Plots.mm)

    # Add legend
    legend_plot = create_legend(
        strategy_pair, color_map, colors, cmap, results["winner_grids"], D_subset;
        legend_fontsize=s_legend_fs, rect_scale=rect_scale,
        legend_width=effective_width, legend_height=legend_height
    )

    final_plot = plot(subplot_layout, legend_plot,
                     layout=grid(2, 1, heights=[subplot_frac, 1 - subplot_frac]),
                     size=(effective_width, total_height),
                     dpi=dpi)

    savefig(final_plot, save_path)
    println("Phase diagram saved to: $(save_path)")
    println("Format: $(uppercase(splitext(save_path)[2][2:end]))")
    display(final_plot)

    return final_plot
end

# ==============================================================================
# Filename Generation and Save Functions
# ==============================================================================

"""
Generate intelligently named filename for competition phase diagram plots.
"""
function generate_phase_diagram_filename(strategy_pair, D_values;
                                         format="pdf", output_dir="plot_results")
    # Map strategy pair to readable names
    strategy_pair_display = if strategy_pair == "0-1"
        "N-T"
    elseif strategy_pair == "1-2"
        "T-P"
    elseif strategy_pair == "0-2"
        "N-P"
    else
        strategy_pair
    end

    # Summarize D range from file data
    D_min = @sprintf("%.3f", minimum(D_values))
    D_max = @sprintf("%.3f", maximum(D_values))
    n_D = length(D_values)

    filename = "competition_phase_diagram_$(strategy_pair_display)_D-$(D_min)-$(D_max)_n$(n_D).$(format)"

    # Ensure output directory exists
    if !isdir(output_dir)
        mkpath(output_dir)
        println("Created output directory: $(output_dir)")
    end

    return joinpath(output_dir, filename)
end

"""
Generate and save competition phase diagram plot as PDF.
"""
function save_phase_diagram(results; format="pdf", output_dir="plot_results",
                            max_plots=9, figure_width=1772,
                            subplot_title_fontsize=22, axis_label_fontsize=20,
                            tick_fontsize=18, legend_fontsize=18,
                            show_subplot_labels=true, subplot_label_fontsize=22,
                            title_gap=0, dpi=300, subplot_titles=nothing)
    println("\n" * "="^70)
    println("Saving competition phase diagram plot...")
    println("="^70)

    strategy_pair = results["strategy_pair"]

    # Extract actual D values from file data
    file_D_values = Float64[]
    for fd in results["file_data_list"]
        append!(file_D_values, fd["D_values"])
    end
    file_D_values = sort(unique(file_D_values))
    D_subset = file_D_values[1:min(max_plots, length(file_D_values))]

    save_path = generate_phase_diagram_filename(strategy_pair, D_subset;
                                                 format=format, output_dir=output_dir)

    plot_robustness_multi_subplots(results;
        save_path=save_path,
        max_plots=max_plots,
        figure_width=figure_width,
        subplot_title_fontsize=subplot_title_fontsize,
        axis_label_fontsize=axis_label_fontsize,
        tick_fontsize=tick_fontsize,
        legend_fontsize=legend_fontsize,
        show_subplot_labels=show_subplot_labels,
        subplot_label_fontsize=subplot_label_fontsize,
        title_gap=title_gap,
        dpi=dpi,
        subplot_titles=subplot_titles)

    println("Saved: $(basename(save_path))")
    println("Output directory: $(output_dir)")
    println("="^70)

    return save_path
end

# ==============================================================================
# Main
# ==============================================================================

function main()
    println("=== Phase Diagram Multi-Subplot Analysis ===")

    # ── Specify your CSV files here ──────────────────────────────────
    csv_files = [
    # Add your CSV file paths here, e.g.:
    # "path/to/competition_summary_pair_12_D_0.04_grid9_combined.csv",
    ]


    if isempty(csv_files)
        println("No CSV files specified. Please edit the csv_files list in main().")
        return 1
    end

    try
        results = load_multiple_csv_files(csv_files)

        println("\n=== Creating Phase Diagram ===")

        # Save as PDF with auto-generated filename
        
        custom_titles = nothing
            
            
  
        save_phase_diagram(results; format="pdf", output_dir="plot_results",
                           subplot_titles=
                           custom_titles
                           )

        return 0

    catch e
        println("Error: $(e)")
        println(sprint(showerror, e, catch_backtrace()))
        return 1
    end
end

main()
