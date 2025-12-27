"""
This module provides functionalities for plotting sundials.
"""

import matplotlib.pyplot as plt
import numpy as np

roman_num = {0:"XXIV", 1:"I", 2:"II", 3:"III", 4:"IV", 5:"V", 6:"VI", 7:"VII", 8:"VIII", 9:"IX", 10:"X", 11:"XI", 12:"XII", 13:"XIII", 14:"XIV", 15:"XV", 16:"XVI", 17:"XVII", 18:"XVIII", 19:"XIX", 20:"XX", 21:"XXI", 22: "XXII", 23:"XXIII"}

class Plotter:
    def __init__(self, xsize=10, ysize=10, cm_per_unit=2):
        self.layout = {}
        self.fig = plt.figure(figsize=(xsize*cm_per_unit/2.54, ysize*cm_per_unit/2.54))
        self.ax = self.fig.add_axes((0,0,1,1))
        self.ax.axis(xmin=-xsize/2, xmax=xsize/2, ymin=-ysize/2, ymax=ysize/2)
        self.xsize = xsize
        self.ysize = ysize
        
    def save(self, filename, resolution=300):
        self.fig.savefig(filename, dpi=resolution)
        
    def text(self, text, x, y, rotation=0, fontsize=12, color="black"):
        self.ax.text(x, y, text, ha="center", va="center", rotation=rotation, fontsize=fontsize, color=color)
        
    def plot_point(self, pos, style=None):
        if style is None:
            style = {}
        style.setdefault("ms", 4)
        style.setdefault("color", "black")
        self.ax.plot([pos[0]], [pos[1]], "o", markersize=style["ms"], color=style["color"])
        
    def find_text_position_straight_line(self, A, B, side="left", position=0.8, position_type="relative", tightness=0.25, rotate=False, rotation_type="perpendicular"):
        A = np.asarray(A)
        B = np.asarray(B)
        
        vec = B-A
        length = np.linalg.norm(vec)
        vec /= length
        
        if np.isclose(vec[0], 0):
            nus = np.array([(-self.ysize/2-A[1])/vec[1], (self.ysize/2-A[1])/vec[1]])
        elif np.isclose(vec[1], 0):
            nus = np.array([(-self.xsize/2-A[0])/vec[0], (self.xsize/2-A[0])/vec[0]])
        else:
            temp_nus = np.array([(-self.xsize/2-A[0])/vec[0], (self.xsize/2-A[0])/vec[0], (-self.ysize/2-A[1])/vec[1], (self.ysize/2-A[1])/vec[1]])
            nus = np.array([])
            for i in range(len(temp_nus)):
                if (i < 2 and (-self.ysize/2-1e-8 <= A[1]+vec[1]*temp_nus[i] <= self.ysize/2+1e-8)) or (i >= 2 and (-self.xsize/2-1e-8 <= A[0]+vec[0]*temp_nus[i] <= self.xsize/2+1e-8)):
                    nus = np.append(nus, temp_nus[i])
           
        if len(nus) > 0:
            if 0 < np.min(nus) < length:
                A = A+np.min(nus)*vec
            if 0 < np.max(nus) < length:
                B = A+np.max(nus)*vec
           
        if rotate:
            if rotation_type == "perpendicular":
                angle = (np.arctan2(vec[1], vec[0])-np.pi/2)*180/np.pi
            elif rotation_type == "parallel_up":
                angle = (np.arctan2(vec[1], vec[0]))*180/np.pi
            elif rotation_type == "parallel_down":
                angle = (np.arctan2(vec[1], vec[0]))*180/np.pi+180
            else:
                raise Exception("rotation_type must be 'perpendicular', 'parallel_up' or 'parallel_down'")             
        else:
            angle = 0
            
        if position_type == "relative":
            lambd = np.linalg.norm(B-A)*position
        elif position_type == "absolute":
            lambd = position
        else:
            raise Exception("position_type must be 'relative' or 'absolute'")
            
        if side == "left":
            coord = A+lambd*vec-tightness*np.array([vec[1], -vec[0]])
        elif side == "right":
            coord = A+lambd*vec+tightness*np.array([vec[1], -vec[0]])
        elif side == "center":
            coord = A+lambd*vec
        else:
            raise Exception("side must be 'left', right' or 'center'")
        return [coord, angle]
    
    def find_text_position_curve(self, X, Y, side="left", position=0.8, position_type="relative", tightness=0.25, rotate=True, rotation_type="perpendicular"):
        X = np.asarray(X)
        Y = np.asarray(Y)
        
        ind_notnan = np.where(np.isfinite(X)&np.isfinite(Y))[0]
        X = X[ind_notnan]
        Y = Y[ind_notnan]
        
        ind_inside = np.where(((-self.xsize/2<X)&(X<self.xsize/2))&((-self.ysize/2<Y)&(Y<self.ysize/2)))[0]
        X = X[ind_inside]
        Y = Y[ind_inside]
        
        dist = np.cumsum(np.append(np.array([0]), np.sqrt((np.roll(X, -1)-X)**2+(np.roll(Y, -1)-Y)**2)[0:-1]))
        
        if np.size(X) < 3:
            return [np.array([np.nan, np.nan]), np.nan]
        
        if position_type == "relative":
            ind = max(min(np.argmin(np.abs(dist-dist[-1]*position)), np.size(X)-2), 1)
        elif position_type == "absolute":
            ind = max(min(np.argmin(np.abs(dist-position)), np.size(X)-2), 1)
        else:
            raise Exception("position_type must be 'relative' or 'absolute'")
            
        tang = np.array([X[ind+1]-X[ind-1], Y[ind+1]-Y[ind-1]])
        tang /= np.linalg.norm(tang)
        
        if rotate:
            if rotation_type == "perpendicular":
                angle = (np.arctan2(tang[1], tang[0])-np.pi/2)*180/np.pi
            elif rotation_type == "parallel_up":
                angle = (np.arctan2(tang[1], tang[0]))*180/np.pi
            elif rotation_type == "parallel_down":
                angle = (np.arctan2(tang[1], tang[0]))*180/np.pi+180
            else:
                raise Exception("rotation_type must be 'perpendicular', 'parallel_up' or 'parallel_down'")
        else:
            angle = 0
                
        if side == "left":
            coord = np.array([X[ind], Y[ind]])-tightness*np.array([tang[1], -tang[0]])
        elif side == "right":
            coord = np.array([X[ind], Y[ind]])+tightness*np.array([tang[1], -tang[0]])
        elif side == "center":
            coord = np.array([X[ind], Y[ind]])
        else:
            raise Exception("side must be 'left', right' or 'center'")
            
        return [coord, angle]
    
    def plot_straight_line(self, A, B, style=None, label=None, labelstyle=None):
        A = np.asarray(A)
        B = np.asarray(B)  
        
        if style is None:
            style = {}
        style.setdefault("lw", 2)
        style.setdefault("color", "blue")
        style.setdefault("linestyle", "-")
        
        if labelstyle is None:
            labelstyle = {}
        labelstyle.setdefault("side", "left")
        labelstyle.setdefault("position", 0.8)
        labelstyle.setdefault("position_type", "relative")
        labelstyle.setdefault("tightness", 0.25)
        labelstyle.setdefault("fontsize", 8)
        labelstyle.setdefault("color", style["color"])
        labelstyle.setdefault("rotate", False)
        labelstyle.setdefault("rotation_type", "perpendicular")
        
        self.ax.plot([A[0], B[0]], [A[1], B[1]], lw=style["lw"], color=style["color"], linestyle=style["linestyle"])
        
        if (label is not None) and np.all(np.isfinite(np.append(A, B))):
            pos = self.find_text_position_straight_line(A, B, side=labelstyle["side"], position=labelstyle["position"], position_type=labelstyle["position_type"], tightness=labelstyle["tightness"], rotate = labelstyle["rotate"], rotation_type=labelstyle["rotation_type"])
            if (-self.xsize/2 < pos[0][0] < self.xsize/2) and (-self.ysize/2 < pos[0][1] < self.ysize/2):
                self.text(label, pos[0][0], pos[0][1], rotation=pos[1], fontsize=labelstyle["fontsize"], color=labelstyle["color"])
            
    def plot_curve(self, X, Y, style=None, label=None, labelstyle=None):
        if style is None:
            style = {}
        style.setdefault("lw", 2)
        style.setdefault("color", "black")
        style.setdefault("linestyle", "-")
        
        if labelstyle is None:
            labelstyle = {}
        labelstyle.setdefault("side", "left")
        labelstyle.setdefault("position", 0.8)
        labelstyle.setdefault("position_type", "relative")
        labelstyle.setdefault("tightness", 0.25)
        labelstyle.setdefault("fontsize", 8)
        labelstyle.setdefault("color", style["color"])
        labelstyle.setdefault("rotate", True)
        labelstyle.setdefault("rotation_type", "perpendicular")
        
        self.ax.plot(X, Y, lw=style["lw"], color=style["color"], linestyle=style["linestyle"])
        
        if label is not None:
            pos = self.find_text_position_curve(X, Y, side=labelstyle["side"], position=labelstyle["position"], position_type=labelstyle["position_type"], tightness=labelstyle["tightness"], rotate=labelstyle["rotate"], rotation_type=labelstyle["rotation_type"])
            if np.all(np.isfinite(pos[0])) and np.isfinite(pos[1]):
                self.text(label, pos[0][0], pos[0][1], rotation=pos[1], fontsize=labelstyle["fontsize"], color=labelstyle["color"])
