import numpy as np
from scipy.spatial.transform import Rotation as spr
import matplotlib.pyplot as plt
import matplotlib.collections
#import matplotlib.mplot3D
# import mayavi.mlab
from Modded_FullCrystalDataSet import Fullcrystaldata


#INPUT:   ===============================
Input_Crystal = "CUB" #ONLY BCT2 for this test (c>a)
# of the style "FCC" or "cF" (capitalization matters!)
dim = (0.8,0.5,01.0,0.6,0.4,0.9)
# print(dim)
# The dimentions (a,b,c,alpha,beta,gamma). Only those defined in a particular crystal structure will be used be the program. Angles are in Radians.

# Crystal definitions taken from https://arxiv.org/pdf/1004.2974.pdf


# Individual Crystal Data ====================================
#imported from another program (used soley to save space with the crystal id library)
def t3D_kplot_builder(dim,Input_Crystal):

    (Corners,Segments,Plot_List,Corner_Labels,b1,b2,b3,Ghosted_Segments,Outer_Segments,unfold_instructions) = Fullcrystaldata(dim,Input_Crystal)

    # Ghosted_Segments = [[1,3],[2,6],[4,7]]
    # Outer_Segments = [[2,6],[2,7],[3,6],[4,7]]
    # unfold_instructions = [["mirror",[0,2,5]],["mirror",[0,3,6]],["mirror",[0,3,1]],["2fmirror",[[0,5,2],[0,3,4]]]]

    Plot_Lines = []
    for seg in Segments:
        Plot_Lines.append([Corners[seg[0]],Corners[seg[1]]])

    Ghosted_Lines = []
    for seg in Ghosted_Segments:
        Ghosted_Lines.append([Corners[seg[0]],Corners[seg[1]]])

    Outer_Lines = []
    for seg in Outer_Segments:
        Outer_Lines.append([Corners[seg[0]],Corners[seg[1]]])



    ###### Rotation:
        
    # rot1 = spr.from_euler('z', 90, degrees=True)
    # Plot_Mirror = []
    # for seg in Plot_Lines:
    #     [m,n] = seg
    #     mrot = np.array(spr.apply(rot1,m))
    #     nrot = np.array(spr.apply(rot1,n))
    #     Plot_Mirror.append([mrot,nrot])
    # Plot_Lines.extend(Plot_Mirror)

    ##### Reflection:
    def mirroronplane(Seg_List,mirror_points):
        [c1,c2,c3] = mirror_points
        p1 = Corners[c1]
        p2 = Corners[c2]
        p3 = Corners[c3]
        plane = [1,2,1,4]
        crossprod = np.cross(p3-p1,p2-p1)
        plane[0:3] = crossprod
        plane[3] = np.dot(crossprod, p3)
        squaredsum = plane[0]**2+plane[1]**2+plane[2]**2
        xyzvec = np.array(plane[0:3])
        normlen = squaredsum**-0.5
        normvec = xyzvec*normlen
        opoint = xyzvec*plane[3]/squaredsum
        def mirrororor(point):
            R = point - np.dot((point-opoint),normvec)*normvec
            mirrorpoint = point + 2*(R-point)
            return mirrorpoint
        Plot_Mirror = []
        for seg in Seg_List:
            [m,n] = seg
            mrot = np.array(mirrororor(m))
            nrot = np.array(mirrororor(n))
            Plot_Mirror.append([mrot,nrot])
        return Plot_Mirror

    def mirroredmirroronplane(Seg_List,sub_mirror,mirror_points):
        [c1,c2,c3] = sub_mirror
        p1 = Corners[c1]
        p2 = Corners[c2]
        p3 = Corners[c3]
        plane = [1,2,1,4]
        crossprod = np.cross(p3-p1,p2-p1)
        plane[0:3] = crossprod
        plane[3] = np.dot(crossprod, p3)
        squaredsum = plane[0]**2+plane[1]**2+plane[2]**2
        xyzvec = np.array(plane[0:3])
        normlen = squaredsum**-0.5
        normvec = xyzvec*normlen
        opoint = xyzvec*plane[3]/squaredsum
        def mirrororor(point):
            R = point - np.dot((point-opoint),normvec)*normvec
            mirrorpoint = point + 2*(R-point)
            return mirrorpoint
        [c4,c5,c6] = mirror_points
        p1a = mirrororor(Corners[c4])
        p2a = mirrororor(Corners[c5])
        p3a = mirrororor(Corners[c6])
        plane = [1,2,1,4]
        crossprod = np.cross(p3a-p1a,p2a-p1a)
        plane[0:3] = crossprod
        plane[3] = np.dot(crossprod, p3a)
        squaredsum = plane[0]**2+plane[1]**2+plane[2]**2
        xyzvec = np.array(plane[0:3])
        normlen = squaredsum**-0.5
        normvec = xyzvec*normlen
        opoint = xyzvec*plane[3]/squaredsum
        def mirrororor(point):
            R = point - np.dot((point-opoint),normvec)*normvec
            mirrorpoint = point + 2*(R-point)
            return mirrorpoint
        Plot_Mirror = []
        for seg in Seg_List:
            [m,n] = seg
            mrot = np.array(mirrororor(m))
            nrot = np.array(mirrororor(n))
            Plot_Mirror.append([mrot,nrot])
        return Plot_Mirror


    for n in unfold_instructions:
        if n[0] == "mirror":
            Outer_Lines.extend(mirroronplane(Outer_Lines,n[1]))
        if n[0] == "2fmirror":
            Outer_Lines.extend(mirroredmirroronplane(Outer_Lines,n[1][0],n[1][1]))


    Single_Outer_Lines = []                 #eliminates redundant outer lines from the mirroring
    for n in Outer_Lines:
        Add = True
        for m in Single_Outer_Lines:
            if all(n[0] == m[0]):
                if all(n[1] == m[1]):
                    Add = False
        if Add == True:
            Single_Outer_Lines.append(n)
    
    return (Single_Outer_Lines,Plot_Lines,Ghosted_Lines,Corners,Corner_Labels)

Single_Outer_Lines,Plot_Lines,Ghosted_Lines,Corners,Corner_Labels = t3D_kplot_builder(dim,Input_Crystal)

ax3d = plt.subplot(1,1,1,projection='3d')

for n in range(np.shape(Single_Outer_Lines)[0]):
    ax3d.plot([Single_Outer_Lines[n][0][0],Single_Outer_Lines[n][1][0]], \
              [Single_Outer_Lines[n][0][1],Single_Outer_Lines[n][1][1]],
              [Single_Outer_Lines[n][0][2],Single_Outer_Lines[n][1][2]], color = 'black', linewidth = 0.8)

for n in range(np.shape(Plot_Lines)[0]):
    ax3d.plot([Plot_Lines[n][0][0],Plot_Lines[n][1][0]],[Plot_Lines[n][0][1],Plot_Lines[n][1][1]],\
              [Plot_Lines[n][0][2],Plot_Lines[n][1][2]],color = 'red')

for n in range(np.shape(Ghosted_Lines)[0]):
    ax3d.plot([Ghosted_Lines[n][0][0],Ghosted_Lines[n][1][0]],[Ghosted_Lines[n][0][1],Ghosted_Lines[n][1][1]],\
              [Ghosted_Lines[n][0][2],Ghosted_Lines[n][1][2]],linestyle='dashed', color = 'red')

#Corner Dots and labels
ax3d.scatter(Corners[:,0],Corners[:,1],Corners[:,2],color = 'red')
for n in range(len(Corner_Labels)):
    ax3d.text(Corners[n][0],Corners[n][1],Corners[n][2],Corner_Labels[n])

#squaring the plot box
[x0,x1] = ax3d.axes.get_xbound()
[y0,y1] = ax3d.axes.get_ybound()
[z0,z1] = ax3d.axes.get_zbound()
xaxismag = np.abs(x1-x0)
yaxismag = np.abs(y1-y0)
zaxismag = np.abs(z1-z0)
largestaxismag = np.max([xaxismag,yaxismag,zaxismag])
ax3d.axes.set_xbound(x0+xaxismag/2-largestaxismag/2,x1-xaxismag/2+largestaxismag/2)
ax3d.axes.set_ybound(y0+yaxismag/2-largestaxismag/2,y1-yaxismag/2+largestaxismag/2)
ax3d.axes.set_zbound(z0+zaxismag/2-largestaxismag/2,z1-zaxismag/2+largestaxismag/2)


plt.show()
    
# exit()   #----------------------------------------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# #Plotting it in Mayavi=============================================

# mayavi.mlab.points3d(Points[:,0],Points[:,1],Points[:,2])
# # mayavi.mlab.points3d(Points[6,0],Points[6,1],Points[6,2],color=(1,0,0))
# # mayavi.mlab.points3d(DOSPoints[:,0],DOSPoints[:,1],DOSPoints[:,2],color=(0,1,0))
# mayavi.mlab.show()
    


