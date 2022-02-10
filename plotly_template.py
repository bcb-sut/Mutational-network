import plotly.graph_objects as go
from plotly.subplots import make_subplots


def make_multiple_plots(plots, rows, cols, types, titles, xlabels, ylabels):
    # specs = [{'type': t} for t in types]
    specs = []
    cnt = 0
    for i in range(rows):
        spec_row = []
        for col in range(cols):
            spec_row.append({'type': types[cnt]})
            cnt += 1
        specs.append(spec_row)

    fig = make_subplots(rows=rows, cols=cols, specs=specs,
                        vertical_spacing=0.05, horizontal_spacing=0.1,
                        subplot_titles=titles)

    idx = 0
    for row in range(1, rows + 1):
        for col in range(1, cols + 1):
            fig.add_trace(plots[idx], row=row, col=col)
            fig.update_yaxes(title_text=ylabels[idx], row=row, col=col)
            fig.update_xaxes(title_text=xlabels[idx], row=row, col=col)
            idx += 1
    return fig


def scatter(x, y):
    return go.Scatter(x=x, y=y,
                      mode='markers',
                      name='original')


def table(header, cells):
    table_ = go.Table(
        header=dict(values=header,
                    fill_color='paleturquoise',
                    align='center'),
        cells=dict(values=cells,
                   fill_color='lavender',
                   align='center'))

    fig = go.Figure(data=[table_])
    return fig


def histogram(x):
    return go.Histogram(x=x, name=None)


def show_go(fig):
    go.Figure(fig).show()


def pie(vals, labels, colors=None):
    pie = go.Pie(values=vals, labels=labels)
    fig = go.Figure(data=[pie])
    if colors:
        fig.update_traces(marker=dict(colors=colors))
    return fig
