#!/usr/bin/env python3
"""
ASCII Art Plotting Library - Matplotlib-style plots using ASCII characters
"""

import math
from typing import List, Tuple, Optional, Union


class ASCIIPlot:
    """Create matplotlib-style plots using ASCII art."""
    
    def __init__(self, width: int = 80, height: int = 24):
        """
        Initialize ASCII plot.
        
        Args:
            width: Width of the plot in characters
            height: Height of the plot in characters
        """
        self.width = width
        self.height = height
        self.title = ""
        self.xlabel = ""
        self.ylabel = ""
        self.data_series = []  # List of (x, y, marker, label) tuples
        self.show_grid = True
        
    def plot(self, x: List[float], y: List[float], marker: str = '•', label: str = ''):
        """
        Add a line/scatter plot.
        
        Args:
            x: X coordinates
            y: Y coordinates
            marker: Character to use for plotting (e.g., '*', 'o', '•', '+')
            label: Label for legend
        """
        if len(x) != len(y):
            raise ValueError("x and y must have the same length")
        self.data_series.append((x, y, marker, label))
        
    def scatter(self, x: List[float], y: List[float], marker: str = 'o', label: str = ''):
        """
        Add a scatter plot.
        
        Args:
            x: X coordinates
            y: Y coordinates
            marker: Character to use for plotting
            label: Label for legend
        """
        self.plot(x, y, marker, label)
        
    def set_title(self, title: str):
        """Set plot title."""
        self.title = title
        
    def set_xlabel(self, label: str):
        """Set x-axis label."""
        self.xlabel = label
        
    def set_ylabel(self, label: str):
        """Set y-axis label."""
        self.ylabel = label
        
    def grid(self, show: bool = True):
        """Toggle grid display."""
        self.show_grid = show
        
    def _get_data_bounds(self) -> Tuple[float, float, float, float]:
        """Get min/max bounds for all data series."""
        if not self.data_series:
            return 0, 1, 0, 1
            
        all_x = []
        all_y = []
        for x, y, _, _ in self.data_series:
            all_x.extend(x)
            all_y.extend(y)
            
        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)
        
        # Add padding
        x_range = x_max - x_min if x_max != x_min else 1
        y_range = y_max - y_min if y_max != y_min else 1
        
        x_min -= x_range * 0.05
        x_max += x_range * 0.05
        y_min -= y_range * 0.05
        y_max += y_range * 0.05
        
        return x_min, x_max, y_min, y_max
    
    def _format_number(self, num: float, width: int = 8) -> str:
        """Format a number for axis labels."""
        if abs(num) < 0.01 or abs(num) >= 10000:
            s = f"{num:.2e}"
        elif abs(num) < 1:
            s = f"{num:.4f}"
        elif abs(num) < 100:
            s = f"{num:.2f}"
        else:
            s = f"{num:.1f}"
        return s.rjust(width)
    
    def _create_canvas(self) -> List[List[str]]:
        """Create an empty canvas."""
        return [[' ' for _ in range(self.width)] for _ in range(self.height)]
    
    def _draw_axes_and_grid(self, canvas: List[List[str]], 
                           x_min: float, x_max: float, 
                           y_min: float, y_max: float,
                           plot_left: int, plot_right: int,
                           plot_top: int, plot_bottom: int):
        """Draw axes, grid, and labels."""
        # Draw axes
        for i in range(plot_top, plot_bottom + 1):
            canvas[i][plot_left] = '│'
        for j in range(plot_left, plot_right + 1):
            canvas[plot_bottom][j] = '─'
        canvas[plot_bottom][plot_left] = '└'
        
        # Draw grid
        if self.show_grid:
            num_grid_lines = 4
            for k in range(1, num_grid_lines):
                # Vertical grid lines
                x_pos = plot_left + k * (plot_right - plot_left) // num_grid_lines
                for i in range(plot_top, plot_bottom):
                    if canvas[i][x_pos] == ' ':
                        canvas[i][x_pos] = '┊'
                
                # Horizontal grid lines
                y_pos = plot_top + k * (plot_bottom - plot_top) // num_grid_lines
                for j in range(plot_left + 1, plot_right + 1):
                    if canvas[y_pos][j] == ' ' or canvas[y_pos][j] == '┊':
                        canvas[y_pos][j] = '┈'
        
        # Y-axis labels
        num_y_ticks = 5
        for k in range(num_y_ticks):
            y_canvas = plot_bottom - k * (plot_bottom - plot_top) // (num_y_ticks - 1)
            y_value = y_min + k * (y_max - y_min) / (num_y_ticks - 1)
            label = self._format_number(y_value, width=8)
            # Place label to the left of axis
            for idx, char in enumerate(label):
                if plot_left - len(label) + idx >= 0:
                    canvas[y_canvas][plot_left - len(label) + idx] = char
        
        # X-axis labels
        num_x_ticks = 5
        for k in range(num_x_ticks):
            x_canvas = plot_left + k * (plot_right - plot_left) // (num_x_ticks - 1)
            x_value = x_min + k * (x_max - x_min) / (num_x_ticks - 1)
            label = self._format_number(x_value, width=8).strip()
            # Center label under tick
            start_pos = x_canvas - len(label) // 2
            for idx, char in enumerate(label):
                if 0 <= start_pos + idx < self.width and plot_bottom + 1 < self.height:
                    canvas[plot_bottom + 1][start_pos + idx] = char
    
    def _plot_data(self, canvas: List[List[str]],
                   x_min: float, x_max: float,
                   y_min: float, y_max: float,
                   plot_left: int, plot_right: int,
                   plot_top: int, plot_bottom: int):
        """Plot data points on the canvas."""
        plot_width = plot_right - plot_left
        plot_height = plot_bottom - plot_top
        
        markers = ['•', '*', '+', 'x', 'o', '◆', '■', '▲']
        
        for series_idx, (x, y, marker, label) in enumerate(self.data_series):
            if not marker:
                marker = markers[series_idx % len(markers)]
            
            for xi, yi in zip(x, y):
                # Convert data coordinates to canvas coordinates
                x_canvas = plot_left + int((xi - x_min) / (x_max - x_min) * plot_width)
                y_canvas = plot_bottom - int((yi - y_min) / (y_max - y_min) * plot_height)
                
                # Check bounds
                if plot_left <= x_canvas <= plot_right and plot_top <= y_canvas <= plot_bottom:
                    canvas[y_canvas][x_canvas] = marker
    
    def show(self) -> str:
        """Generate and return the ASCII plot as a string."""
        if not self.data_series:
            return "No data to plot"
        
        # Get data bounds
        x_min, x_max, y_min, y_max = self._get_data_bounds()
        
        # Calculate plot area (leaving room for labels)
        plot_left = 10  # Space for y-axis labels
        plot_right = self.width - 2
        plot_top = 2 if self.title else 1
        plot_bottom = self.height - 3  # Space for x-axis labels
        
        # Create canvas
        canvas = self._create_canvas()
        
        # Draw axes and grid
        self._draw_axes_and_grid(canvas, x_min, x_max, y_min, y_max,
                                plot_left, plot_right, plot_top, plot_bottom)
        
        # Plot data
        self._plot_data(canvas, x_min, x_max, y_min, y_max,
                       plot_left, plot_right, plot_top, plot_bottom)
        
        # Add title
        if self.title:
            title_start = (self.width - len(self.title)) // 2
            for idx, char in enumerate(self.title):
                if 0 <= title_start + idx < self.width:
                    canvas[0][title_start + idx] = char
        
        # Add axis labels
        if self.ylabel:
            # Vertical text for y-label
            y_label_col = 0
            y_label_row = (plot_top + plot_bottom) // 2 - len(self.ylabel) // 2
            for idx, char in enumerate(self.ylabel):
                if 0 <= y_label_row + idx < self.height:
                    canvas[y_label_row + idx][y_label_col] = char
        
        if self.xlabel:
            x_label_row = self.height - 1
            x_label_start = (self.width - len(self.xlabel)) // 2
            for idx, char in enumerate(self.xlabel):
                if 0 <= x_label_start + idx < self.width:
                    canvas[x_label_row][x_label_start + idx] = char
        
        # Add legend
        legend_labels = [label for _, _, _, label in self.data_series if label]
        if legend_labels:
            legend_row = plot_top
            legend_col = plot_right - 20
            for idx, (_, _, marker, label) in enumerate(self.data_series):
                if label and legend_col > plot_left:
                    legend_text = f"{marker} {label}"
                    for char_idx, char in enumerate(legend_text):
                        if 0 <= legend_col + char_idx < self.width:
                            canvas[legend_row + idx][legend_col + char_idx] = char
        
        # Convert canvas to string
        return '\n'.join(''.join(row) for row in canvas)
    
    def display(self):
        """Print the plot to console."""
        print(self.show())


def demo():
    """Demonstrate the ASCII plotting capabilities."""
    import math
    
    # Example 1: Sine and Cosine
    print("Example 1: Sine and Cosine waves")
    print("=" * 80)
    x = [i * 0.1 for i in range(63)]
    y_sin = [math.sin(xi) for xi in x]
    y_cos = [math.cos(xi) for xi in x]
    
    plot1 = ASCIIPlot(width=80, height=24)
    plot1.set_title("Trigonometric Functions")
    plot1.set_xlabel("x (radians)")
    plot1.set_ylabel("y")
    plot1.plot(x, y_sin, marker='*', label='sin(x)')
    plot1.plot(x, y_cos, marker='+', label='cos(x)')
    plot1.display()
    
    print("\n" * 2)
    
    # Example 2: Quadratic function
    print("Example 2: Quadratic Function")
    print("=" * 80)
    x2 = [i * 0.2 - 5 for i in range(51)]
    y2 = [xi**2 - 2*xi - 3 for xi in x2]
    
    plot2 = ASCIIPlot(width=80, height=20)
    plot2.set_title("y = x² - 2x - 3")
    plot2.set_xlabel("x")
    plot2.set_ylabel("y")
    plot2.plot(x2, y2, marker='•')
    plot2.display()
    
    print("\n" * 2)
    
    # Example 3: Scatter plot
    print("Example 3: Scatter Plot with Multiple Series")
    print("=" * 80)
    x3a = [1, 2, 3, 4, 5, 6, 7, 8]
    y3a = [2, 4, 3, 5, 7, 6, 8, 9]
    x3b = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]
    y3b = [3, 5, 4, 6, 5, 7, 8]
    
    plot3 = ASCIIPlot(width=80, height=20)
    plot3.set_title("Data Comparison")
    plot3.set_xlabel("Time (s)")
    plot3.set_ylabel("Value")
    plot3.scatter(x3a, y3a, marker='◆', label='Series A')
    plot3.scatter(x3b, y3b, marker='■', label='Series B')
    plot3.display()
    
    print("\n" * 2)
    
    # Example 4: Exponential
    print("Example 4: Exponential Growth")
    print("=" * 80)
    x4 = [i * 0.1 for i in range(31)]
    y4 = [math.exp(xi * 0.5) for xi in x4]
    
    plot4 = ASCIIPlot(width=80, height=20)
    plot4.set_title("Exponential: y = e^(0.5x)")
    plot4.set_xlabel("x")
    plot4.set_ylabel("y")
    plot4.plot(x4, y4, marker='*')
    plot4.display()


if __name__ == "__main__":
    demo()
