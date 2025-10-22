import pandas as pd
import numpy as np
import plotly
import seaborn
#import plotly.offline as pyo
import plotly.graph_objs as go


def offset_signal(signal, marker_offset):
    if abs(signal) <= marker_offset:
        return 0
    return signal - marker_offset if signal > 0 else signal + marker_offset

def lollipop_plot(dataframe_genes):

    # Offset the line length by the marker size to avoid overlapping
    marker_offset = 0.04


    data = [
        go.Scatter(
            x=list(range(len(dataframe_genes['global_score']))),
            y=list(dataframe_genes['global_score'].values),
            text=dataframe_genes['name'],
            mode='markers',
            marker=dict(color='red')
        )
    ]

    # Use the 'shapes' attribute from the layout to draw the vertical lines
    layout = go.Layout(
        shapes=[dict(
            type='line',
            xref='x',
            yref='y',
            x0=i,
            y0=0,
            x1=i,
            y1=offset_signal(list(dataframe_genes['global_score'].values)[i], marker_offset)
            line=dict(
                color='grey',
                width=1
            )
        ) for i in range(len(dataframe_genes['global_score']))],
        title='Predicted association scores', xaxis_title='Index', yaxis_title='Score'
    )

    # Plot the chart


    fig = go.Figure(data, layout)
    fig.update_layout(yaxis_range=[0,1], xaxis_range=[-1,50.5])
    pyo.iplot(fig)


def plot_model_loss(x_test, y_test, x_train, y_train):
    layout = go.Layout(title='Model loss per epoch', xaxis_title='Epoch', yaxis_title='Loss')
    data = [
        go.Scatter(
            x=x_train,
            y=y_train,
            mode='lines',name='Training set'),
        go.Scatter(
            x=x_test,
            y=y_test,
            mode='lines',name='Test set'),            
    ]
    fig = go.Figure(data,layout)
    pyo.iplot(fig)
    return

def plot_accuracy_loss(x_test, y_test, x_train, y_train):
    layout = go.Layout(title='Model accuracy per epoch', xaxis_title='Epoch', yaxis_title='Loss')
    data = [
        go.Scatter(
            x=x_train,
            y=y_train,
            mode='lines',name='Training set'),
        go.Scatter(
            x=x_test,
            y=y_test,
            mode='lines',name='Test set'),            
    ]
    fig = go.Figure(data,layout)
    pyo.iplot(fig)
    return

def plot_matrix(true_targets, predicted_targets, matrix):
    """
    The matrix should be a numpy array in this shape
    array([[2, 2, 2, ..., 2, 0, 1],
       [3, 1, 1, ..., 0, 1, 0],
       [0, 3, 1, ..., 1, 0, 1],
       ...,
       [0, 1, 1, ..., 0, 1, 2],
       [1, 2, 3, ..., 1, 0, 1],
       [1, 1, 0, ..., 1, 0, 0]], shape=(7, 180))

    """
    data=[go.Heatmap(z=matrix,x=predicted_targets, y=true_targets,colorscale='Viridis')]
    layout = go.Layout(title='Confidence in target prediction', xaxis_title='True target', yaxis_title='Predicted target')
    fig = go.Figure(data,layout)
    pyo.iplot(fig)   
    return
