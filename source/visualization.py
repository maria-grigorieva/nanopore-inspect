import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd

def plot_distribution(sequences, smoothing, mode='proportion'):
    fig = go.Figure()
    fig_title = ("Positional distribution of occurrences proportional to the total number "
                 "of sequences (Proportion) and Reads") if mode == 'proportion' else \
        ("Distribution of positional occurrences of sequences")

    for p in sequences:
        title = p['type']
        if 'value_counts' in p:
            df = pd.DataFrame(p['value_counts'])
            total_proportion = p['total_proportion']
            total_reads = p['total_reads']
            noise = p['noise_level']
            fig_legend = (f'{title}: Total Proportion = {total_proportion}, '
                          f'number of reads = {total_reads}, level of noise = {noise}') if mode == 'proportion' else \
                (f'{title}: Total Reads = {total_reads}, proportion = {total_proportion}, level of noise = {noise}')

            fig.add_trace(go.Scatter(
                x=df['index'],
                y=df[mode],
                mode='lines',
                name=fig_legend,
                yaxis='y1'  # Specify the first Y-axis
                )
            )

            if len(sequences) == 1:
                # Add dots at specific x positions (for proportions)
                peak_positions = [int(i['peak_index']) for i in p['peaks']]
                fig.add_trace(go.Scatter(
                    x=peak_positions,
                    y=[df.loc[df['index'] == x, mode].values[0] for x in peak_positions],
                    text=[i[f'total_{mode}'] for i in p['peaks']],
                    textposition='top center',
                    textfont=dict(color='red'),
                    mode='markers+text',
                    marker=dict(size=8, color='red'),
                    showlegend=False
                ))

                # Smoothing line
                if smoothing != 'None':
                    x_smooth = df['index']
                    y_smooth = df['smoothed']
                    fig.add_trace(go.Scatter(
                        x=x_smooth,
                        y=y_smooth,
                        mode='lines',
                        showlegend=False,
                        line=dict(color='red')
                    ))

    # Update layout with dual Y-axis (proportion on left, reads on right)
    fig.update_layout(
        width=1200,
        height=600,
        legend=dict(
            x=0,
            y=1.0,
            xanchor='left',
            yanchor='bottom'
        ),
        barmode='overlay',
        title={
            'text': fig_title,
            'y': 0.05,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'bottom'
        },
        # Define the primary Y-axis for Proportion
        yaxis=dict(
            title=mode,
            side="left",
            showgrid=True
        )
    )
    # Set opacity for the traces
    fig.update_traces(opacity=0.75)

    return fig

def make_peaks_subplots(sequences):
    specs = [[{"type": "table"}] for i in range(0, len(sequences))]
    titles = []
    for item in sequences:
        seq_type = item['type']
        if 'average_peaks_distancev' in item:
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