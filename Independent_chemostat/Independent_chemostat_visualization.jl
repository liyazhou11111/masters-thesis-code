using CSV, DataFrames, Plots, Printf, Statistics
using Plots.PlotMeasures
using ColorSchemes

# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

function auto_detect_intervals(parameter_values)
    """Automatically detects the interval between points in a parameter range."""
    if length(parameter_values) < 2; return 1.0; end
    return parameter_values[2] - parameter_values[1]
end

function smart_tick_selection(values; max_ticks=8)
    """Selects a subset of values to use as ticks to avoid overcrowding an axis."""
    if length(values) <= max_ticks; return values; end
    step = max(1, div(length(values), max_ticks))
    return values[1:step:end]
end

function get_subplot_label(idx)
    """Generate subplot label (a), (b), (c), etc."""
    return "(" * string(Char('a' + idx - 1)) * ")"
end

function add_subplot_label!(p, label; fontsize=28, x_offset=-0.25, y_offset=1.1)
    """Add a subplot label outside the axis area (above-left corner)."""
    xl = xlims(p)
    yl = ylims(p)
    x_pos = xl[1] + x_offset * (xl[2] - xl[1])
    y_pos = yl[1] + y_offset * (yl[2] - yl[1])
    annotate!(p, [(x_pos, y_pos, text(label, :left, :bottom, fontsize, :bold))])
end

function create_independent_legend(strategy_pair, active_categories;
                                   legend_fontsize=18, rect_scale=0.6, legend_width=1772, legend_height=200)
    """Creates a custom legend plot showing only categories present in the data."""

    # Define the color palette for the legend
    colors = [
        RGB(0.5, 0.5, 0.5),      # 1: Extinct - dark gray
        RGB(0.95, 0.8, 0.25),    # 2: Strategy N - orange/yellow
        RGB(0.35, 0.8, 0.9),     # 3: Strategy T - blue
        RGB(0.95, 0.6, 0.5),     # 4: Strategy P - red/coral
        RGB(0.4, 0.75, 0.55),    # 5: Stable equal - green
        RGB(0.6, 0.7, 0.75),     # 6: Neutral equal - gray-blue
        RGB(0.7, 0.5, 0.8)       # 7: Multi-stability - purple
    ]

    # Build descriptive labels based on strategy pair
    name1, name2 = if strategy_pair == "0-1"
        ("N", "T")
    elseif strategy_pair == "1-2"
        ("T", "P")
    else
        ("S" * split(strategy_pair, "-")[1], "S" * split(strategy_pair, "-")[2])
    end

    all_categories = Dict(
        1 => "Extinct",
        2 => "$name1 > $name2",
        3 => if strategy_pair == "0-1" "$name2 > $name1" else "$name1 > $name2" end,
        4 => "$name2 > $name1",
        5 => "$name1 ≈ $name2 (stable)",
        6 => "$name1 ≈ $name2 (neutral)",
        7 => "Multi-stability"
    )

    # Filter to only categories present in the data
    categories = Dict(k => v for (k, v) in all_categories if k in active_categories)

    # Create an empty plot canvas for the legend
    legend_plot = plot(xlims=(0,1), ylims=(0,1), framestyle=:none, grid=false, legend=false, size=(legend_width, legend_height))

    n_items = length(categories)
    if n_items == 0; return legend_plot; end

    # Calculate positions for each legend item
    x_positions = range(0.08, 0.92, length=n_items)
    categories_list = sort(collect(categories), by=x->x[1])

    # Draw each legend item (a colored box and a label)
    for (idx, (color_idx, label)) in enumerate(categories_list)
        x_center = x_positions[idx]
        exact_color = colors[color_idx]
        rect_width = 0.8 / n_items * rect_scale
        rect_height = 0.2

        plot!(legend_plot,
              Shape([x_center - rect_width/2, x_center + rect_width/2, x_center + rect_width/2, x_center - rect_width/2],
                    [0.6 - rect_height/2, 0.6 - rect_height/2, 0.6 + rect_height/2, 0.6 + rect_height/2]),
              color=exact_color, linecolor=:black, linewidth=1)

        annotate!(legend_plot, x_center, 0.25, text(label, :center, legend_fontsize))
    end

    return legend_plot
end

# ==============================================================================
#  CORE DATA AND PLOTTING FUNCTIONS
# ==============================================================================

function load_and_organize_data(csv_files, column_map)
    """
    Loads data from multiple CSV files, combines them, and organizes the data into grids
    for plotting, using a specific column map.
    """
    println("=== Loading Independent Chemostat Grid Data ===")
    all_dfs = DataFrame[]
    
    # Handle single file input
    if csv_files isa String
        csv_files = [csv_files]
    end
    
    for (i, csv_file) in enumerate(csv_files)
        if !isfile(csv_file)
            println("⚠️ Warning: Skipping missing file: $(csv_file)")
            continue
        end
        println("Loading file $(i)/$(length(csv_files)): $(basename(csv_file))")
        
        try
            # Force the 'strategy_pair' and 'outcome' columns to be read as text
            df = CSV.read(csv_file, DataFrame; 
                          stringtype=String, 
                          dateformat=nothing,
                          types=Dict(column_map[:strategy_pair] => String,
                                    column_map[:outcome] => String))

            # Check if all required columns exist in this file
            required_cols = values(column_map)
            if !all(c -> c in names(df), required_cols)
                println("  - ❌ ERROR: File $(basename(csv_file)) is missing required columns. Skipping.")
                println("    Required: $(collect(required_cols))")
                println("    Found:    $(names(df))")
                continue
            end
            push!(all_dfs, df)
        catch e
            println("  - ❌ ERROR: Could not read file $(basename(csv_file)): $(e). Skipping.")
        end
    end
    if isempty(all_dfs)
        error("No valid CSV files were loaded.")
    end
    
    # Combine all loaded dataframes into one
    combined_df = vcat(all_dfs...)
    println("Successfully loaded and combined $(length(all_dfs)) files.")
    
    # --- Organize and clean the combined data ---
    println("\n--- Organizing data... ---")

    # Rename user-specified columns to a standard internal format
    for (std_name, user_name) in column_map
        if user_name in names(combined_df)
            rename!(combined_df, Symbol(user_name) => std_name)
        end
    end

    # Clean the 'strategy_pair' column by replacing date-like strings (legacy fix)
    if any(s -> occursin(r"\d{4}-\d{2}-\d{2}", s), combined_df.strategy_pair)
        println("Fixing date-like strings in 'strategy_pair' column...")
        combined_df.strategy_pair = replace.(combined_df.strategy_pair, r"\d{4}-01-01" => "0-1")
        combined_df.strategy_pair = replace.(combined_df.strategy_pair, r"\d{4}-02-01" => "1-2")
    end
    println("Strategy pairs found: $(unique(combined_df.strategy_pair))")
    println("Outcomes found: $(unique(combined_df.outcome))")

    # Extract parameter ranges
    T1_values = sort(unique(combined_df.T1))
    T2_values = sort(unique(combined_df.T2))
    D_values = sort(unique(combined_df.D))
    strategy_pairs = sort(unique(combined_df.strategy_pair))

    # Define outcome to color index mapping
    # Color indices: 1=Extinct, 2=N, 3=T, 4=P, 5=Stable_equal, 6=Neutral_equal, 7=Multi-stability
    outcome_to_color = Dict(
        "Both_extinct" => 1,
        "Strategy_0_higher_biomass" => 2,
        "Strategy_1_higher_biomass" => 3,
        "Strategy_2_higher_biomass" => 4,
        "Stable_equal" => 5,
        "Neutral_equal" => 6,
        "Multi-stability" => 7,
        "No_Data" => 1  # Treat as extinct/missing
    )

    # Create data grids for each D value and strategy pair
    outcome_grids = Dict()
    for D in D_values, pair in strategy_pairs
        key = "$(D)_$(pair)"
        outcome_grid = fill(1, length(T1_values), length(T2_values))  # Default to 1 (Extinct/Unknown)
        
        subset = combined_df[(combined_df.D .== D) .& (combined_df.strategy_pair .== pair), :]
        for row in eachrow(subset)
            t1_idx = findfirst(x -> x == row.T1, T1_values)
            t2_idx = findfirst(x -> x == row.T2, T2_values)
            if t1_idx !== nothing && t2_idx !== nothing
                outcome_grid[t1_idx, t2_idx] = get(outcome_to_color, row.outcome, 1)
            end
        end
        outcome_grids[key] = outcome_grid
    end
    println("--- Data organization complete ---")

    return Dict(
        "T1_values" => T1_values, 
        "T2_values" => T2_values, 
        "D_values" => D_values, 
        "strategy_pairs" => strategy_pairs, 
        "outcome_grids" => outcome_grids
    )
end

function plot_independent_multi_subplots(results;
                                       strategy_pair="1-2",
                                       save_path="independent_outcome.png",
                                       max_plots=9,
                                       # 150mm @ 300dpi ≈ 1772 pixels (matching LaTeX textwidth)
                                       figure_width=1772,
                                       # Font size parameters
                                       subplot_title_fontsize=22,
                                       label_fontsize=20,
                                       tick_fontsize=18,
                                       legend_fontsize=18,
                                       # Subplot label parameters
                                       show_subplot_labels=true,
                                       subplot_label_fontsize=22,
                                       title_gap=0,
                                       dpi=300)
    """
    Creates and saves a multi-panel heatmap from the organized grid data.
    Includes flexible layouts, a global legend, configurable font sizes,
    and optional subplot labels (a), (b), (c), etc.
    """
    # --- Setup ---
    T1_values = results["T1_values"]
    T2_values = results["T2_values"]
    D_values = results["D_values"]
    D_subset = D_values[1:min(max_plots, length(D_values))]
    n_plots = length(D_subset)

    # Determine subplot layout automatically based on the number of plots
    if n_plots == 1
        layout_tuple = (1, 1)
    elseif n_plots == 2
        layout_tuple = (1, 2)
    elseif n_plots <= 3
        layout_tuple = (1, 3)
    elseif n_plots <= 4
        layout_tuple = (2, 2)
    elseif n_plots <= 6
        layout_tuple = (2, 3)
    elseif n_plots <= 8
        layout_tuple = (2, 4)
    else
        layout_tuple = (3, 3)
    end

    # Calculate figure height based on layout (fixed A4 width, proportional height)
    n_rows, n_cols = layout_tuple
    row_height = round(Int, figure_width / n_cols * 0.85)
    figure_height = n_rows * row_height
    legend_height = 200

    # Generate title gap string (creates space between title and plot for subplot labels)
    gap_str = repeat("\n", title_gap)

    println("\n=== Creating Plot: $(save_path) ===")
    println("Strategy pair: $(strategy_pair), Number of subplots: $(n_plots)")
    
    # --- Axes and Color Setup ---
    selected_T2_values = smart_tick_selection(T2_values; max_ticks=4)
    selected_T1_values = smart_tick_selection(T1_values; max_ticks=4)
    T1_boundaries = [T1_values[1] - 0.5*auto_detect_intervals(T1_values); [v + 0.5*auto_detect_intervals(T1_values) for v in T1_values]]
    T2_boundaries = [T2_values[1] - 0.5*auto_detect_intervals(T2_values); [v + 0.5*auto_detect_intervals(T2_values) for v in T2_values]]

    # Define colors: Extinct, N, T, P, Stable_equal, Neutral_equal, Multi-stability
    colors = [
        RGB(0.5, 0.5, 0.5),      # 1: Extinct - dark gray
        RGB(0.95, 0.8, 0.25),    # 2: Strategy N - orange/yellow
        RGB(0.35, 0.8, 0.9),     # 3: Strategy T - blue
        RGB(0.95, 0.6, 0.5),     # 4: Strategy P - red/coral
        RGB(0.4, 0.75, 0.55),    # 5: Stable equal - green
        RGB(0.6, 0.7, 0.75),     # 6: Neutral equal - gray-blue
        RGB(0.75, 0.55, 0.2)       # 7: Multi-stability - purple
    ]

    # --- Subplot Generation ---
    subplots = []
    for (idx, D) in enumerate(D_subset)
        key = "$(D)_$(strategy_pair)"
        heat_data = get(results["outcome_grids"], key, fill(1, length(T1_values), length(T2_values)))

        # Determine if this subplot should show axis labels
        show_xlabel = (idx-1) >= (layout_tuple[1]-1) * layout_tuple[2]
        show_ylabel = (idx-1) % layout_tuple[2] == 0

        # Create base plot with axes and styling (no heatmap interpolation)
        p = plot(
            title=@sprintf("D = %.3f", D) * gap_str,
            xlabel=show_xlabel ? "T₂ (S₁,ᵢₙ OFF Duration)" : "",
            ylabel=show_ylabel ? "T₁ (S₁,ᵢₙ ON Duration)" : "",
            framestyle=:box,
            xticks=(selected_T2_values, [isinteger(v) ? string(Int(v)) : string(v) for v in selected_T2_values]),
            yticks=selected_T1_values,
            xlims=(T2_boundaries[1], T2_boundaries[end]),
            ylims=(T1_boundaries[1], T1_boundaries[end]),
            legend=false,
            titlefontsize=subplot_title_fontsize,
            xlabelfontsize=label_fontsize,
            ylabelfontsize=label_fontsize,
            xtickfontsize=tick_fontsize,
            ytickfontsize=tick_fontsize
        )

        # Draw each cell as a solid-colored rectangle (no interpolation)
        n_T1, n_T2 = size(heat_data)
        for j in 1:n_T2
            for i in 1:n_T1
                c = colors[clamp(heat_data[i, j], 1, length(colors))]
                plot!(p, Shape(
                    [T2_boundaries[j], T2_boundaries[j+1], T2_boundaries[j+1], T2_boundaries[j]],
                    [T1_boundaries[i], T1_boundaries[i], T1_boundaries[i+1], T1_boundaries[i+1]]
                ), fillcolor=c, linecolor=c, linewidth=0, label=false)
            end
        end

        # Add subplot label (a), (b), (c), etc. outside the axis area (top-left corner)
        if show_subplot_labels
            add_subplot_label!(p, get_subplot_label(idx), fontsize=subplot_label_fontsize)
        end

        push!(subplots, p)
    end
    
    # Fill remaining slots with empty plots if needed
    while length(subplots) < prod(layout_tuple)
        push!(subplots, plot(framestyle=:none))
    end
    
    # --- Final Plot Assembly ---
    subplot_layout = plot(subplots...; layout=layout_tuple,
                          left_margin=12Plots.mm, right_margin=8Plots.mm,
                          bottom_margin=10Plots.mm, top_margin=8Plots.mm)

    # Collect all category indices that actually appear in the data
    active_categories = Set{Int}()
    for D in D_subset
        key = "$(D)_$(strategy_pair)"
        grid = get(results["outcome_grids"], key, nothing)
        if grid !== nothing
            union!(active_categories, unique(grid))
        end
    end

    # Add legend at the bottom (only showing categories present in the data)
    legend_plot = create_independent_legend(strategy_pair, active_categories; legend_fontsize=legend_fontsize, legend_width=figure_width, legend_height=legend_height)
    final_plot = plot(subplot_layout, legend_plot,
                      layout=@layout([A{0.88h}; B{0.12h}]),
                      size=(figure_width, figure_height + legend_height),
                      dpi=dpi)
    
    savefig(final_plot, save_path)
    println("✅ Plot saved to: $(save_path)")
    display(final_plot)
    
    return final_plot
end

# ==============================================================================
#  FILENAME GENERATION AND SAVE FUNCTIONS
# ==============================================================================

function generate_phase_diagram_filename(strategy_pair, D_values;
                                         format="pdf", output_dir="plot_results")
    """
    Generate intelligently named filename for phase diagram plots.

    Parameters:
    - strategy_pair: Strategy pair string (e.g., "0-1", "1-2")
    - D_values: Array of D values included in the plot
    - format: File format (default "pdf")
    - output_dir: Output directory (default "plot_results")

    Returns:
    - Full file path
    """
    # Map strategy pair to readable names
    strategy_names = if strategy_pair == "0-1"
        "N-T"
    elseif strategy_pair == "1-2"
        "T-P"
    else
        "S$(split(strategy_pair, "-")[1])-S$(split(strategy_pair, "-")[2])"
    end

    # Summarize D range
    D_min = @sprintf("%.3f", minimum(D_values))
    D_max = @sprintf("%.3f", maximum(D_values))
    n_D = length(D_values)

    filename = "independent_phase_diagram_$(strategy_names)_D-$(D_min)-$(D_max)_n$(n_D).$(format)"

    # Ensure output directory exists
    if !isdir(output_dir)
        mkpath(output_dir)
        println("📁 Created output directory: $(output_dir)")
    end

    return joinpath(output_dir, filename)
end

function save_phase_diagrams(results; format="pdf", output_dir="plot_results",
                             max_plots=9, figure_width=1772,
                             subplot_title_fontsize=22, label_fontsize=20,
                             tick_fontsize=18, legend_fontsize=18,
                             show_subplot_labels=true, subplot_label_fontsize=22,
                             title_gap=0, dpi=300)
    """
    Generate and save all phase diagram plots as PDF files.

    Parameters:
    - results: Data dictionary from load_and_organize_data
    - format: File format (default "pdf")
    - output_dir: Output directory (default "plot_results")
    - Other parameters: Passed through to plot_independent_multi_subplots
    """
    println("\n" * "="^70)
    println("🎨 Saving all independent phase diagram plots...")
    println("="^70)

    saved_files = String[]
    D_values = results["D_values"]
    D_subset = D_values[1:min(max_plots, length(D_values))]

    for (idx, pair) in enumerate(results["strategy_pairs"])
        println("\n📊 [$(idx)/$(length(results["strategy_pairs"]))] Saving phase diagram: strategy pair $(pair)...")
        try
            save_path = generate_phase_diagram_filename(pair, D_subset;
                                                        format=format, output_dir=output_dir)
            plot_independent_multi_subplots(results;
                strategy_pair=pair,
                save_path=save_path,
                max_plots=max_plots,
                figure_width=figure_width,
                subplot_title_fontsize=subplot_title_fontsize,
                label_fontsize=label_fontsize,
                tick_fontsize=tick_fontsize,
                legend_fontsize=legend_fontsize,
                show_subplot_labels=show_subplot_labels,
                subplot_label_fontsize=subplot_label_fontsize,
                title_gap=title_gap,
                dpi=dpi)
            push!(saved_files, save_path)
            println("✅ Saved: $(basename(save_path))")
        catch e
            println("❌ Failed to save phase diagram for pair $(pair): $e")
        end
    end

    # Print summary
    println("\n" * "="^70)
    println("🎉 All phase diagrams saved successfully!")
    println("="^70)
    println("📁 Output directory: $(output_dir)")
    println("📄 File format: $(uppercase(format))")
    println("📊 Number of files saved: $(length(saved_files))")
    println("\nSaved files:")
    for (i, file) in enumerate(saved_files)
        println("  $i. $(basename(file))")
    end
    println("="^70)

    return saved_files
end

# ==============================================================================
#  MAIN EXECUTION BLOCK
# ==============================================================================

function main()
    println("=== Independent Chemostat Grid Outcome Visualization ===")

    # --------------------------------------------------------------------------
    #  ACTION REQUIRED 1: Configure your column names here
    # --------------------------------------------------------------------------
    column_map = Dict(
        :T1 => "T1",
        :T2 => "T2", 
        :D => "D",
        :strategy_pair => "strategy_pair",
        :outcome => "outcome",
        :winner => "winner",
        :consensus => "consensus_value"
    )

    # --------------------------------------------------------------------------
    #  ACTION REQUIRED 2: List your CSV files here
    # --------------------------------------------------------------------------
    # Example for single file:
    # csv_files = "/path/to/your/file.csv"
    
    # Example for multiple files:
     csv_files = [
    # Add your CSV file paths here, e.g.:
    # "path/to/independent_grid_combined_summary_D_0.02_0-1.csv",
     ]


    
    if isempty(csv_files)
        println("⚠️ Please add your CSV file paths to the `csv_files` variable in the `main` function.")
        return
    end
    
    try
        # Load and process data
        results = load_and_organize_data(csv_files, column_map)
        
        # Generate and save plots as PDF for each strategy pair
        save_phase_diagrams(results; format="pdf", output_dir="plot_results")
        
    catch e
        println("\n❌ An error occurred during the process.")
        println(sprint(showerror, e, catch_backtrace()))
        return 1
    end
    return 0
end

# Run the main function
main()