import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd

def plot_distribution(sequences, smoothing, mode='proportion'):
    # Define a modern color palette (colorblind-friendly)
    colors = [
        '#2E86AB',  # Blue
        '#A23B72',  # Purple
        '#3C8C51',  # Green
        '#F18F01',  # Orange
        '#C73E1D',  # Red
        '#7B4B94'  # Deep Purple
    ]

    # Define consistent font settings
    TITLE_FONT_SIZE = 24
    AXIS_FONT_SIZE = 16
    LEGEND_FONT_SIZE = 14
    PEAK_LABEL_FONT_SIZE = 12
    FONT_FAMILY = "Arial, Helvetica, sans-serif"

    fig = go.Figure()

    # Determine figure title based on mode
    fig_title = ("Positional Distribution of Occurrences<br>"
                 "<sup>Proportional to Total Number of Sequences (Proportion) and Reads</sup>") if mode == 'proportion' else \
        ("Distribution of Positional Occurrences<br>"
         "<sup>Sequence Analysis Results</sup>")

    for idx, p in enumerate(sequences):
        color = colors[idx % len(colors)]
        title = p['type']

        if 'value_counts' in p:
            df = pd.DataFrame(p['value_counts'])
            total_proportion = f"{p['total_proportion']:.3f}"
            total_reads = f"{p['total_reads']:,}"
            noise = f"{p['noise_level']:.3f}"

            # Create more readable legend labels
            # fig_legend = (f'{title}<br>'
            #               f'Total Prop: {total_proportion} | '
            #               f'Reads: {total_reads}<br>'
            #               f'Noise: {noise}') if mode == 'proportion' else \
            #     (f'{title}<br>'
            #      f'Reads: {total_reads} | '
            #      f'Prop: {total_proportion}<br>'
            #      f'Noise: {noise}')
            fig_legend = (f'{title}:  '
                          f'proportion {total_proportion} | reads {total_reads} ') \
                if mode == 'proportion' else \
                (f'{title}:  '
                 f'{total_reads} | {total_proportion}')

            # Main curve with enhanced styling
            fig.add_trace(go.Scatter(
                x=df['index'],
                y=df[mode],
                mode='lines',
                name=fig_legend,
                line=dict(
                    color=color,
                    width=3
                ),
                yaxis='y1'
            ))

            # Peak values with improved visibility
            peak_positions = [int(i['peak_index']) for i in p['peaks']]
            peak_values = [df.loc[df['index'] == x, mode].values[0] for x in peak_positions]
            peak_labels = [f"{i[f'total_{mode}']:.3f}" for i in p['peaks']]

            fig.add_trace(go.Scatter(
                x=peak_positions,
                y=peak_values,
                text=peak_labels,
                textposition='top center',
                textfont=dict(
                    size=PEAK_LABEL_FONT_SIZE,
                    color=color,
                    family=FONT_FAMILY
                ),
                mode='markers+text',
                marker=dict(
                    size=10,
                    color=color,
                    symbol='diamond',
                    line=dict(
                        color='white',
                        width=1
                    )
                ),
                showlegend=False
            ))

            # Smoothed line with refined styling
            if smoothing != 'None':
                fig.add_trace(go.Scatter(
                    x=df['index'],
                    y=df['smoothed'],
                    mode='lines',
                    line=dict(
                        color=color,
                        dash='dot',
                        width=2
                    ),
                    opacity=0.5,
                    showlegend=False
                ))

    # Enhanced layout with modern styling
    fig.update_layout(
        # Size and margins
        width=1200,
        height=700,
        margin=dict(
            l=80,
            r=50,
            t=100,
            b=80
        ),

        # Plot style
        plot_bgcolor='white',
        paper_bgcolor='white',

        # Title
        title=dict(
            text=fig_title,
            y=0.95,
            x=0.5,
            xanchor='center',
            yanchor='top',
            font=dict(
                size=TITLE_FONT_SIZE,
                family=FONT_FAMILY
            )
        ),

        # Legend
        legend=dict(
            x=0,
            y=1.0,
            xanchor='left',
            yanchor='bottom',
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='lightgray',
            borderwidth=1,
            font=dict(
                size=LEGEND_FONT_SIZE,
                family=FONT_FAMILY
            )
        ),

        # Axes
        xaxis=dict(
            title=dict(
                text="Position",
                font=dict(
                    size=AXIS_FONT_SIZE,
                    family=FONT_FAMILY
                )
            ),
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='lightgray',
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True
        ),

        yaxis=dict(
            title=dict(
                text=mode.capitalize(),
                font=dict(
                    size=AXIS_FONT_SIZE,
                    family=FONT_FAMILY
                )
            ),
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='lightgray',
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True
        ),

        # Hover behavior
        # hovermode='x unified'
    )

    # Add custom hover template
    fig.update_traces(
        hovertemplate="<b>Position:</b> %{x}<br>" +
                      f"<b>{mode.capitalize()}:</b> %{{y:.3f}}<br>" +
                      "<extra></extra>"
    )

    return fig


def make_peaks_subplots(sequences):
    specs = [[{"type": "table"}] for i in range(0, len(sequences))]
    titles = []
    for item in sequences:
        seq_type = item['type']
        if 'average_peaks_distance' in item:
            dist = item['average_peaks_distance']
            titles.append(f'{seq_type}: average distance between peaks is {round(dist)}')
        else:
            titles.append(seq_type)
    fig = make_subplots(
        rows=len(sequences), cols=1,
        shared_xaxes=True,
        vertical_spacing=0.03,
        specs=specs,
        subplot_titles=titles
        # subplot_titles=[i['type'] for i in sequences]
    )

    for i, item in enumerate(sequences, start=1):
        sequence_type = item['type']
        peaks_df = pd.DataFrame(item['peaks'])
        selected_columns = ['peak_index', 'left_bases', 'right_bases', 'total_proportion', 'total_occurrences',
                            'peak_dist']
        # Get the existing columns in peaks_df that are also in selected_columns
        existing_columns = [col for col in peaks_df.columns if col in selected_columns]

        peaks_df = peaks_df[existing_columns]

        fig.add_trace(
            go.Table(
            header=dict(values=peaks_df.columns),
            cells=dict(values=peaks_df.transpose().values)
        ),
            row=i, col=1
        )

    fig.update_layout(width=1000, height=1000, title_text='Subplots')
    return fig