
from __future__ import annotations
import functools
from itertools import chain

import tkinter as tk
from tkinter import ttk, SE
from tkinter import messagebox
from tkinter.filedialog import asksaveasfilename
import os

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
import matplotlib

import numpy as np

from json_interface import load_lattice
import central_equation
import potentials
from density import plot_densities

# tell matplotlib to use TkAgg so we can show plots on tkinter windows
matplotlib.use("TkAgg")


#####################################################################################################
# define a tooltip widget for use throughout the gui

class tooltip(ttk.Label):
    """Tkinter widget for adding hover-over tooltips to other widgets."""
    def __init__(self, target: ttk.Widget, **kwargs):
        ttk.Label.__init__(self, root, **kwargs, wraplength = 100)
        self.target = target
        target.bind('<Enter>', self.show)
        target.bind('<Leave>', self.on_leave_target)
        self.bind('<Leave>', self.on_leave_self)

    def show(self, e: tk.Event):
        # enter event coordinates are relative to the widget being entered
        # thus, we add its "global" position to compensate
        # winfo_rootx and y give screen space coordinates, so must subtract window position
        # winfo_x and y give coordinates relative to parent
        offset_x = self.target.winfo_rootx() - root.winfo_rootx()
        offset_y = self.target.winfo_rooty() - root.winfo_rooty()
        self.place(x = e.x + offset_x, y = e.y + offset_y, anchor = SE)

    def on_leave_self(self, e: tk.Event):
        x, y = root.winfo_pointerxy()
        widget_under_mouse = root.winfo_containing(x, y)
        if widget_under_mouse is not self.target:
            self.hide()

    def on_leave_target(self, e: tk.Event):
        x, y = root.winfo_pointerxy()
        widget_under_mouse = root.winfo_containing(x, y)
        if widget_under_mouse is not self:
            self.hide()

    def hide(self):
        self.place_forget()

#####################################################################################################
# set up whole window

root = tk.Tk()
root.title("Lattice Bands")
root.geometry("1000x600")

pop = tk.Toplevel(root)
pop.destroy()

"""Tkinker Popup Addition RF ================================================================================================"""
def lattice_popup():                #Builds popup
    global pop
    if pop.winfo_exists() == True:
        pop.destroy()
    pop = tk.Toplevel(root)
    pop.title("Lattice Geometries")
    # pop.geometry("900x400+600+200")
    fig = Figure(figsize = (7, 9))  

    scrollbar = tk.Scrollbar(pop)
    scrollbar.grid(column=1, sticky='nw')

    lattice_index_that_one_number = np.where(CrystalNames==lattice.get())[0]            #Getting crystal plotting data
    Input_Crystal = CrystalList[lattice_index_that_one_number,1]
    dave = CrystalList[lattice_index_that_one_number,2:8][0].tolist()
    dim = (float(dave[0]),float(dave[1]),float(dave[2]),float(dave[3]),float(dave[4]),float(dave[5]))
    modnumber = CrystalList[lattice_index_that_one_number,8]

    ax3d2 = fig.add_subplot(211,projection='3d')        #Real Space Plot

    from Modded_FullCrystalDataSet import Fullcrystaldata              #Jimmied in alternate data set from the Json method
    (Corners,Segments,Plot_List,Corner_Labels,b1,b2,b3,Ghosted_Segments,Outer_Segments,unfold_instructions,a1,a2,a3) = Fullcrystaldata(dim,Input_Crystal,modnumber)
    
    Point_Locations = [[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,0],[1,1,1]]
    Arrow_List = [[0,1],[0,2],[0,3]]
    Ghosted_List = [[1,5],[1,6],[2,4],[2,6],[3,4],[3,5],[4,7],[5,7],[6,7]]

    Real_Corners = []                       #This whole section is building the real space plot
    for n in Point_Locations:
        Real_Corners.append(np.array(a1)*n[0]+np.array(a2)*n[1]+np.array(a3)*n[2])

    Label_Name_YN = True
    msn = 400
    markersizeslist = [msn,msn,msn,msn,msn,msn,msn,msn]
    for n in Real_Corners:
        if Label_Name_YN == True:
            ax3d2.scatter(n[0],n[1],n[2],color = (0.3,0.4,1),sizes=markersizeslist,label = "Atom Positions")
        else:
            ax3d2.scatter(n[0],n[1],n[2],color = (0.3,0.4,1),sizes=markersizeslist)
        Label_Name_YN = False

    Label_Name_YN = True
    for n in Ghosted_List:
        m1 = Real_Corners[n[0]]
        m2 = Real_Corners[n[1]]
        if Label_Name_YN == True:
            ax3d2.plot([m1[0],m2[0]],[m1[1],m2[1]],[m1[2],m2[2]], color = 'black',linestyle='dashed', linewidth = 0.8,label="Primative Unit Cell")
        else:
            ax3d2.plot([m1[0],m2[0]],[m1[1],m2[1]],[m1[2],m2[2]], color = 'black',linestyle='dashed', linewidth = 0.8)
        Label_Name_YN = False
    
    Label_Name_YN = True    
    for n in Arrow_List:
        m1 = Real_Corners[n[0]]
        m2 = Real_Corners[n[1]]   #Arrow goes from m1 to m2
        if Label_Name_YN == True:
            ax3d2.quiver(m1[0],m1[1],m1[2],m2[0]-m1[0],m2[1]-m1[1],m2[2]-m1[2],color=(0,1,1),label="Lattice Vectors")
        else:
            ax3d2.quiver(m1[0],m1[1],m1[2],m2[0]-m1[0],m2[1]-m1[1],m2[2]-m1[2],color=(0,1,1))
        Label_Name_YN = False

    a1pos = Real_Corners[1]-0.5*(Real_Corners[1]-Real_Corners[0])
    a2pos = Real_Corners[2]-0.5*(Real_Corners[2]-Real_Corners[0])
    a3pos = Real_Corners[3]-0.5*(Real_Corners[3]-Real_Corners[0])
    ax3d2.text(a1pos[0],a1pos[1],a1pos[2],"a$_1$")
    ax3d2.text(a2pos[0],a2pos[1],a2pos[2],"a$_2$")
    ax3d2.text(a3pos[0],a3pos[1],a3pos[2],"a$_3$")

    #squaring the plot box
    zscalefactor = 1.4
    [x0,x1] = ax3d2.axes.get_xbound()
    [y0,y1] = ax3d2.axes.get_ybound()
    [z0,z1] = ax3d2.axes.get_zbound()
    xaxismag = np.abs(x1-x0)
    yaxismag = np.abs(y1-y0)
    zaxismag = np.abs(z1-z0)*zscalefactor
    largestaxismag = np.max([xaxismag,yaxismag,zaxismag])
    ax3d2.axes.set_xbound(x0+xaxismag/2-largestaxismag/2,x1-xaxismag/2+largestaxismag/2)
    ax3d2.axes.set_ybound(y0+yaxismag/2-largestaxismag/2,y1-yaxismag/2+largestaxismag/2)
    ax3d2.axes.set_zbound(z0+zaxismag/2/zscalefactor-largestaxismag/2/zscalefactor,\
                        z1-zaxismag/2/zscalefactor+largestaxismag/2/zscalefactor)
    
    ax3d2.legend(loc = 'lower right')
    ax3d2.set_title('Real Space: units: length')
    

    ax3d = fig.add_subplot(212,projection='3d')            #Momentum Space Plot

    from t3D_frame_construction import t3D_kplot_builder            #This uses a bunch of the infatructure from the "testing" plot program

    #Part copied direct from t3D_frame_construction:
    Single_Outer_Lines,Plot_Lines,Ghosted_Lines,Corners,Corner_Labels = t3D_kplot_builder(dim,Input_Crystal,modnumber)

    Label_Name_YN = True
    for n in range(np.shape(Single_Outer_Lines)[0]):
        if Label_Name_YN == True:
            ax3d.plot([Single_Outer_Lines[n][0][0],Single_Outer_Lines[n][1][0]], \
                    [Single_Outer_Lines[n][0][1],Single_Outer_Lines[n][1][1]],
                    [Single_Outer_Lines[n][0][2],Single_Outer_Lines[n][1][2]], color = 'black', linewidth = 0.8,label="Full Brillouin Zone")
        else:
            ax3d.plot([Single_Outer_Lines[n][0][0],Single_Outer_Lines[n][1][0]], \
                    [Single_Outer_Lines[n][0][1],Single_Outer_Lines[n][1][1]],
                    [Single_Outer_Lines[n][0][2],Single_Outer_Lines[n][1][2]], color = 'black', linewidth = 0.8)
        Label_Name_YN = False

    Label_Name_YN = True
    for n in range(np.shape(Plot_Lines)[0]):
        if Label_Name_YN == True:
            ax3d.plot([Plot_Lines[n][0][0],Plot_Lines[n][1][0]],[Plot_Lines[n][0][1],Plot_Lines[n][1][1]],\
                    [Plot_Lines[n][0][2],Plot_Lines[n][1][2]],color = 'red',label="Reduced Brillouin Zone")
        else:
            ax3d.plot([Plot_Lines[n][0][0],Plot_Lines[n][1][0]],[Plot_Lines[n][0][1],Plot_Lines[n][1][1]],\
                    [Plot_Lines[n][0][2],Plot_Lines[n][1][2]],color = 'red')
        Label_Name_YN = False

    Label_Name_YN = True
    for n in range(np.shape(Ghosted_Lines)[0]):
        if Label_Name_YN == True:
            ax3d.plot([Ghosted_Lines[n][0][0],Ghosted_Lines[n][1][0]],[Ghosted_Lines[n][0][1],Ghosted_Lines[n][1][1]],\
                    [Ghosted_Lines[n][0][2],Ghosted_Lines[n][1][2]],linestyle='dashed', color = 'red', label="Non-Sampled Symmetry Points")
        else:
            ax3d.plot([Ghosted_Lines[n][0][0],Ghosted_Lines[n][1][0]],[Ghosted_Lines[n][0][1],Ghosted_Lines[n][1][1]],\
                    [Ghosted_Lines[n][0][2],Ghosted_Lines[n][1][2]],linestyle='dashed', color = 'red')
        Label_Name_YN = False

    #Corner Dots and labels
    ax3d.scatter(Corners[:,0],Corners[:,1],Corners[:,2],color = 'red')
    for n in range(len(Corner_Labels)):
        ax3d.text(Corners[n][0],Corners[n][1],Corners[n][2],Corner_Labels[n])

    #squaring the plot box
    zscalefactor = 1.4
    [x0,x1] = ax3d.axes.get_xbound()
    [y0,y1] = ax3d.axes.get_ybound()
    [z0,z1] = ax3d.axes.get_zbound()
    xaxismag = np.abs(x1-x0)
    yaxismag = np.abs(y1-y0)
    zaxismag = np.abs(z1-z0)*zscalefactor
    largestaxismag = np.max([xaxismag,yaxismag,zaxismag])
    ax3d.axes.set_xbound(x0+xaxismag/2-largestaxismag/2,x1-xaxismag/2+largestaxismag/2)
    ax3d.axes.set_ybound(y0+yaxismag/2-largestaxismag/2,y1-yaxismag/2+largestaxismag/2)
    ax3d.axes.set_zbound(z0+zaxismag/2/zscalefactor-largestaxismag/2/zscalefactor,\
                        z1-zaxismag/2/zscalefactor+largestaxismag/2/zscalefactor)

    #End Copied part.
    
    ax3d.legend(loc="lower right")
    ax3d.set_title("Momentum Space: units: $\pi$/length")

    canvas = FigureCanvasTkAgg(fig, master = pop)   
    canvas.draw() 
    canvas.get_tk_widget().grid(sticky="nesw",column=0,row=0)
    scrollbar["command"] = canvas.get_tk_widget().yview
    # scrollbar.set(0.5,1)
    
#End added popup section

mainframe = ttk.Frame(root, padding = "12")
mainframe.place(relx = 0, rely = 0, relheight = 1, relwidth = 1)
root.columnconfigure(0, weight = 2)
# root.columnconfigure(1, weight = 1)
root.rowconfigure(0, weight = 1)

interfaceframe = ttk.Frame(mainframe, padding = 10)
interfaceframe.grid(column = 1, row = 0, sticky = "NESW")
# interfaceframe.rowconfigure(0, weight = 1)
# interfaceframe.rowconfigure(1, weight = 1)
interfaceframe.rowconfigure(5, weight = 2)





#####################################################################################################
# lattice selector

lattice_label = ttk.Label(interfaceframe, text = "Select lattice to use:")
lattice_label.grid(column = 0, row = 0, sticky = "S", columnspan = 2)

# files: list[str] = os.listdir("lattices")
# json_files = []
# OtherData_names = []

# # TODO: add a name field to the json and load from that?
# for filename in files:
#     if filename.endswith(".json"):
#         json_files.append(filename[:-5])


from Modded_FullCrystalDataSet import CrystalList           #Modded by RF to feed crystal data from my data set instead of the json files
CrystalNames = CrystalList[:,0]

lattice = tk.StringVar()
lattice.set("Cubic (CUB,cP)") #formerly "2d_square"
file_menu = ttk.OptionMenu(interfaceframe, lattice, lattice.get(), *CrystalNames) #formerly json_files
file_menu.grid(column = 0, row = 1, sticky = "N", pady = (0, 20), columnspan = 2)




#####################################################################################################
# define number entry boxes

# fails for negative values
def _validate_int(text: str):
    return text.isdigit() or text == ""

class IntString(tk.StringVar):
    """
    Wrapper for tk.StringVar for use with int_entry.
    Allows entry as a string but ensures and integer value.
    """
    def __init__(self, *args, default: int = 0, **kwargs):
        tk.StringVar.__init__(self, *args, **kwargs)
        self._default = default

    def get(self) -> int:
        tempval = super().get()
        if tempval == "":
            return self._default
        else:
            return int(tempval)

class int_entry(ttk.Entry):
    def __init__(self, *args, **kwargs):
        ttk.Entry.__init__(self, *args, **kwargs)
        validator = self.register(_validate_int)
        # this is the incantaion I found to make text validation work
        self.config(validate = "key", validatecommand = (validator, '%P'))


def _validate_float(text: str):
    # try to split at a decimal point
    match text.split("."):
        # if there isn't a decimal point, it should be an integer
        case [integer]:
            return _validate_int(integer)
        # if there is, both before and after the decimal point should "look like" integers
        case [integer, fraction]:
            return _validate_int(integer) and _validate_int(fraction)
        # otherwise (more than one '.'), fail
        case _:
            return False
            
class FloatString(tk.StringVar):
    """
    Wrapper for tk.StringVar for use with float_entry.
    Allows entry of a string but ensures a float value.
    """
    def __init__(self, *args, default: float = 1, **kwargs):
        tk.StringVar.__init__(self, *args, **kwargs)
        self._default = default

    def get(self) -> float:
        tempval = super().get()
        if tempval == "":
            return self._default
        else:
            return float(tempval)

class float_entry(ttk.Entry):
    def __init__(self, *args, **kwargs):
        ttk.Entry.__init__(self, *args, **kwargs)
        validator = self.register(_validate_float)
        # this is the incantaion I found to make text validation work
        self.config(validate = "key", validatecommand = (validator, '%P'))

#####################################################################################################
# potential selector

potential_frame = ttk.Frame(interfaceframe)
potential_frame.grid(column = 0, row = 5, sticky = "N", columnspan = 2)


strength_label = ttk.Label(interfaceframe, text = "Potential strength:")
strength_label.grid(column = 0, row = 2, sticky = "E")

strength = FloatString(value = "1")
strength_entry = float_entry(interfaceframe, textvariable = strength, width = 4)
strength_entry.grid(column = 1, row = 2, sticky = "W", padx = (10, 0), pady = (0, 10))


potentials_to_plot: list[potential_entry] = []

def add_potential():
    new_pot = potential_entry(potential_frame)
    new_pot.pack(pady = (0, 10))
    potentials_to_plot.append(new_pot)

add_potential_buton = ttk.Button(interfaceframe, text = "Add potential to plot", command = add_potential)
add_potential_buton.grid(column = 0, row = 3, columnspan = 2, sticky = "N", pady = (0, 10))


potential_map = \
{
    "Empty lattice": potentials.empty_v,
    "Simple": potentials.simple_v,
    "Typical": potentials.typical_v,            #Formerlly Called "Quadratic. Changed by Richard"
}
_max_pot_width = max(len(k) for k in potential_map.keys())

line_styles_map = \
{
    "Solid": "-",
    "Dashed": "--",
    "Dotted": ":",
    "Dash-dot": "-."
}
_max_line_width = max(len(k) for k in line_styles_map.keys())

colors = \
[
    "red",
    "blue",
    "green",
    "black",
    "grey",
    "darkred",
    "darkblue",
    "darkgreen",
    "darkgray",
    "rainbow",
]
_max_color_width = max(len(k) for k in colors)

class potential_entry(ttk.Frame):
    """GUI element for displaying and editing parameters of 
    potenitals being plotted, as well as plotting styles.
    """
    # TODO list:
    # Add labels for most elements to explain what they do
    # Maybe check to make sure the delete button doesn't have a memory leak
    def __init__(self, *args, **kwargs):
        ttk.Frame.__init__(self, *args, **kwargs)

        self.separator = ttk.Separator(self, orient = "horizontal")
        self.separator.grid(row = 0, column = 0, columnspan = 100, sticky = "EW", pady = (0, 5))

        self.approximation_label = ttk.Label(self, text = "Select approximation:")
        self.approximation_label.grid(row = 1, column = 0, columnspan = 2, sticky = "W")

        self.potential_name = tk.StringVar(value = "Empty lattice")
        self.potential_selector = ttk.Combobox(self, textvariable = self.potential_name, state = "readonly", values = list(potential_map.keys()), width = _max_pot_width)
        self.potential_selector.current(0)
        self.potential_selector.grid(row = 1, column = 2, columnspan = 3, padx = 5, sticky = "W")

        # self.scale_label = ttk.Label(self, text = "Strength:")
        # self.scale_label.grid(row = 1, column = 2, columnspan = 2, sticky = "E", padx = (5, 0))

        self.scale_var = FloatString(value = "1", default = 1)
        # self.scale_entry = num_entry(self, textvariable = self.scale_var, width = 4)
        # self.scale_entry.grid(row = 1, column = 4, columnspan = 2, sticky = "W", padx = (0, 5))

        self.del_btn = ttk.Button(self, command = self.remove, text = "Remove", width = 8)
        self.del_btn.grid(row = 1, column = 10, sticky = "E")

        self.style_label = ttk.Label(self, text = "Line style:")
        self.style_label.grid(row = 2, column = 0, sticky = "W")

        self.linestyle = tk.StringVar(value = "Solid")
        self.style_selector = ttk.Combobox(self, textvariable = self.linestyle, state = "readonly", values = list(line_styles_map.keys()), width = _max_line_width)
        self.style_selector.grid(row = 3, column = 0, sticky = "W")

        self.color_label = ttk.Label(self, text = "Line color:")
        self.color_label.grid(row = 2, column = 1, sticky = "W", columnspan = 2)

        self.linecolor = tk.StringVar(value = "black")
        self.color_selector = ttk.Combobox(self, textvariable = self.linecolor, state = "readonly", values = colors, width = _max_color_width)
        self.color_selector.grid(row = 3, column = 1, columnspan = 2, sticky = "W")

        self.alpha_label = ttk.Label(self, text = "Line opacity:")
        self.alpha_label.grid(row = 2, column = 4, sticky = "W")

        # use a 0-255 range to leverage existing code for integer handling, but convert to 0-1 later
        self.linealpha = FloatString(value = "1", default = 1)
        self.alpha_entry = float_entry(self, textvariable = self.linealpha, width = 4)
        self.alpha_entry.grid(row = 3, column = 4, sticky = "W")

        self.density_label = ttk.Label(self, text = "Check to make density of states plot:")
        self.density_label.grid(row = 4, column = 0, columnspan = 4)

        self.is_density_checked = tk.BooleanVar(value = False)
        self.density_check = ttk.Checkbutton(self, variable = self.is_density_checked, command = self.density_callback)
        self.density_check.grid(row = 4, column = 10)

    def density_callback(self):
        for entry in potentials_to_plot:
            if entry is not self:
                entry.is_density_checked.set(False)

    def get_potential(self):
        return functools.partial(potential_map[self.potential_name.get()], scale = strength.get())

    def get_line_style(self):
        return line_styles_map[self.linestyle.get()]

    def get_line_color(self):
        return self.linecolor.get()

    def get_line_alpha(self):
        return self.linealpha.get()

    def remove(self):
        self.pack_forget()
        potentials_to_plot.remove(self)

# have something be plotted by default at least
add_potential()
potentials_to_plot[0].is_density_checked.set(True)

#####################################################################################################
# plot parameters - resolution and range

matrix_size_label = ttk.Label(interfaceframe, text = "Size of matrix to use:")
matrix_size_label.grid(column = 0, row = 5, columnspan = 20, sticky = "S")

matrix_size_tip = tooltip(matrix_size_label, text = "The size of the square matrix that will be used to compute the bands. Higher is \"better\" but much slower.")

matrix_size_var = IntString(value = "9", default = 9)
matrix_size_entry = int_entry(interfaceframe, textvariable = matrix_size_var)
matrix_size_entry.grid(column = 0, row = 6, columnspan = 20)

band_count_label = ttk.Label(interfaceframe, text = "Number of bands to plot:")
band_count_label.grid(column = 0, row = 7, columnspan = 20)

band_count_tip = tooltip(band_count_label, text = "The number of bands to plot. The actual number plotted will never be more than the matrix size.")

band_count_var = IntString(value = "9", default = 9)
band_count_entry = int_entry(interfaceframe, textvariable = band_count_var)
band_count_entry.grid(column = 0, row = 8, columnspan = 20)

resolution_label = ttk.Label(interfaceframe, text = "Resolution of plot:", wraplength = 150)
resolution_label.grid(column = 0, row = 9, columnspan = 20)

resolution_tip = tooltip(resolution_label, text = "The number of samples between symmetry points. Higher is smoother but slower.")

resolution_var = IntString(value = "50", default = 50)
resolution_entry = int_entry(interfaceframe, textvariable = resolution_var)
resolution_entry.grid(column = 0, row = 10, columnspan = 20)


#####################################################################################################
# setting up the go button and plotting function
# TODO: look into moving some of these algorithm details to their own function or maybe even module

def plot_bands():
    band_axes.clear()
    density_axes.clear()

    # lattice_filepath = f"lattices/{lattice.get()}.json"
    # lat = load_lattice(lattice_filepath)
    # # vertical_lines = lat.vertical_lines
    # points = lat.points
    # point_names = lat.point_names
    
    lattice_index_that_one_number = np.where(CrystalNames==lattice.get())[0]            
    Input_Crystal = CrystalList[lattice_index_that_one_number,1]
    dave = CrystalList[lattice_index_that_one_number,2:8][0].tolist()
    dim = (float(dave[0]),float(dave[1]),float(dave[2]),float(dave[3]),float(dave[4]),float(dave[5]))
    modnumber = CrystalList[lattice_index_that_one_number,8]
    #Extracting crystal data:
    from Modded_FullCrystalDataSet import Fullcrystaldata              #Jimmied in alternate data set from the Json method
    (Corners,Segments,Plot_List,Corner_Labels,b1,b2,b3,Ghosted_Segments,Outer_Segments,unfold_instructions,a1,a2,a3) = Fullcrystaldata(dim,Input_Crystal,modnumber)
    
    points = []         #Short section to build Ben's varibles out of mine.
    point_names = []
    for n in range(len(Plot_List)):
        if Plot_List[n] != "break":
            points.append(Corners[Plot_List[n]])
            point_names.append(Corner_Labels[Plot_List[n]])

    vertical_lines = range(len(points)) 

    for point in vertical_lines:
        band_axes.axvline(point, linestyle = "--", color = (.5, .5, .5, .5))

    matrix_size = matrix_size_var.get()
    resolution = resolution_var.get()

    band_paths = []
    for n in range(len(points)-1):
        band_paths.append(np.linspace(points[n], points[n + 1], resolution))
    # band_paths = l__at.get_paths(resolution)          #Original method

    plot_ranges = []
    for i in range(len(band_paths)):
        plot_ranges.append(np.linspace(i, i+1, resolution))
    
    path = np.array(list(chain.from_iterable(band_paths)))
    plot_space = np.array(list(chain.from_iterable(plot_ranges)))

    density_energies = None

    for entry in potentials_to_plot:
        potential = entry.get_potential()
    
        style_params = {}
        style_params.update(linestyle = entry.get_line_style())
        if entry.get_line_color() != "rainbow":
            style_params.update(color = to_rgba(entry.get_line_color(), entry.get_line_alpha()))

        energy_bands = central_equation.get_bands(path, potential, matrix_size) [:band_count_var.get()]
        for band in energy_bands:
            band_axes.plot(plot_space, band, **style_params)

        if entry.is_density_checked.get():
            density_energies = energy_bands

    band_axes.set_xlabel("High Symmetry Points", fontsize = "xx-large")
    band_axes.set_xticks(list(range(len(points))))
    band_axes.set_xticklabels(point_names)
    band_axes.set_ylabel(r"Energy, in units of $\frac{ħ^2}{2m}(\frac{π}{a})^2$", fontsize = "xx-large")

    band_axes.relim()
    band_axes.autoscale()

    if density_energies is not None:
        plot_densities(density_energies, density_axes)

        density_axes.set_xlabel("Density", fontsize = "xx-large")
        density_axes.set_xlim(left = 0)
        density_axes.set_visible(True)
    else:
        density_axes.set_visible(False)

    canvas.draw()

go_button = ttk.Button(interfaceframe, text = "Plot bands", command = plot_bands)
go_button.grid(column = 0, row = 20, columnspan = 20, pady = 10)

"""Added Popup Button"""
popup_button = ttk.Button(interfaceframe, text = "Show Lattices", command = lattice_popup)
popup_button.grid(column = 0, row = 21, columnspan = 20, pady = 10)


def print_plot():
    figure_name = asksaveasfilename(initialdir = "figures", defaultextension = ".png")
    canvas.print_png(figure_name)
    
print_button = ttk.Button(interfaceframe, text = "Export Plot as png", command = print_plot)
print_button.grid(column = 0, row = 22, columnspan = 20)

# def print_lattices():
#     figure_name = asksaveasfilename(initialdir = "figures", defaultextension = ".png")
#     lattice_popup.canvas.print_png(figure_name)

# print_button = ttk.Button(interfaceframe, text = "Export Lattices as png", command = print_lattices)
# print_button.grid(column = 0, row = 23, columnspan = 20)


#####################################################################################################
# setting up the matplotlib canvas

# TODO: see if the density subplot can be hidden/resummoned at will

fig = Figure()

band_axes: Axes
density_axes: Axes

# define a gridspec for the displayed subplots, giving more room to the bands than the density plot
grid_spec = \
{
    "width_ratios": (3, 1),
    "wspace": .0,
}
sub_opts = \
{
    "xticklabels": [],
}
band_axes, density_axes = fig.subplots(1, 2, sharey = True, gridspec_kw = grid_spec, subplot_kw = sub_opts)

canvas = FigureCanvasTkAgg(fig, master = mainframe)
canvas.get_tk_widget().grid(column = 0, row = 0, sticky = "NESW")
# make the canvas resize with the window
mainframe.rowconfigure(0, weight = 1)
mainframe.columnconfigure(0, weight = 1)

root.mainloop()
