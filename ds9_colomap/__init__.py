from matplotlib.cm import register_cmap
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from astropy.io import fits as pyfits
from matplotlib.colors import LogNorm
from astropy import wcs
from matplotlib import patches

def draw_reg(ax,file_reg,pixsize):
    f = open(file_reg, 'r')

    start = False
    line = f.readline()
    while line:
        if start == True:
            if 'color' in line:
                cds9 = line.split('color=')[1].split()[0]
                if cds9 == 'white':
                    cc = 'w'
                elif cds9 == 'green':
                    cc = 'lime'
                elif cds9 == 'cyan':
                    cc = 'cyan'
                elif cds9 == 'black':
                    cc = 'k'
                elif cds9 == 'magenta':
                    cc = 'magenta'
                elif cds9 == 'blue':
                    cc = 'b'
                elif cds9 == 'yellow':
                    cc = 'yellow'
                elif cds9 == 'red':
                    cc = 'r'
                else:
                    cc = 'lime'
            else:
                cc = 'lime'
            if 'width' in line:
                ww = float(line.split('width=')[1].split()[0])
            else:
                ww = 1
            temp = line.split('(')
            tipo = temp[0]
            temp = temp[1].split(',')
            xx, yy = float(temp[0]), float(temp[1])
            if f_format == 'fk5':
                pixim = w.all_world2pix([[xx, yy]], 0)
                xx, yy = pixim[0][0], pixim[0][1]
            else:
                xx, yy = xx - 1, yy - 1

            if tipo == 'circle':
                temp2 = temp[2].split(')')[0]
                if f_format == 'fk5':
                    if '\"' in temp2:
                        rr = float(temp2.split('\"')[0]) / pixsize
                    elif '\'' in temp2:
                        rr = float(temp2.split('\'')[0]) * 60 / pixsize
                    else:
                        rr = float(temp2) * 3600 / pixsize
                else:
                    rr = float(temp2)
                circle = plt.Circle((xx, yy), rr, color=cc, fill=False, clip_on=True, lw=ww)
                ax.add_artist(circle)

            elif tipo == 'ellipse':
                tt = float(temp[4].split(')')[0])
                if f_format == 'fk5':
                    temp2 = temp[2]
                    if '\"' in temp2:
                        rx = float(temp2.split('\"')[0]) / pixsize
                    elif '\'' in temp2:
                        rx = float(temp2.split('\'')[0]) * 60 / pixsize
                    else:
                        rx = float(temp2) * 3600 / pixsize
                    temp3 = temp[3]
                    if '\"' in temp3:
                        ry = float(temp3.split('\"')[0]) / pixsize
                    elif '\'' in temp3:
                        ry = float(temp3.split('\'')[0]) * 60 / pixsize
                    else:
                        ry = float(temp3) * 3600 / pixsize
                else:
                    rx = float(temp[2])
                    ry = float(temp[3])
                ellipse = patches.Ellipse((xx, yy), rx, ry, color=cc, angle=tt, linewidth=ww, fill=False)
                ax.add_artist(ellipse)
            else:
                print("ERROR: %s not yet supported, try using circle or ellipse" % tipo)

        if 'fk5' in line:
            start = True
            f_format = 'fk5'
        if 'image' in line:
            start = True
            f_format = 'image'

        line = f.readline()

    f.close()

class ds9_colormap:
    def __init__(self, CIAO='', ds9_cmap=''):
        if CIAO == '':
            if 'ASCDS_CALIB' in os.environ:
                CIAO = os.environ['ASCDS_CALIB'][:-5]  #
                print("Found CIAO: %s in environment variables" % CIAO)
            else:
                raise ValueError("ERROR: CIAO folder not specified: use CIAO=/my/path/to/ciao")
        if not os.path.isdir(CIAO):
            raise ValueError("ERROR: CIAO folder %s does not exists" % CIAO)
        if ds9_cmap == '':
            raise ValueError("ERROR: ds9 colormap not specified: use e.g. ds9_cmap='sls'")

        ds9_cmap_loc = CIAO + "/data/" + ds9_cmap + ".lut"
        if os.path.isfile(ds9_cmap_loc):
            rr, gg, bb = np.loadtxt(ds9_cmap_loc, unpack=True)
            xx = np.linspace(0, 1, len(rr))
            ds9cmap = {'red': lambda v: np.interp(v, xx, rr),
                       'green': lambda v: np.interp(v, xx, gg),
                       'blue': lambda v: np.interp(v, xx, bb)}
            self.ds9_cmap_name = 'ds9%s' % ds9_cmap
            self.ds9_cmap = ds9_cmap
            register_cmap(self.ds9_cmap_name, data=ds9cmap)
            print("ds9 colormap %s registered as 'ds9%s'" % (ds9_cmap, ds9_cmap))
            print(
                        "you can use it directly as: plt.imshow('myimage.fits',cmap='%s',norm=LogNorm(),vmin=7e-7,vmax=1e-3,origin='lower')" % self.ds9_cmap_name)
            print("... to make an image try *.show(image='example.fits',vmin=7e-7,vmax=1e-3,log=True)")
            print("... to display the colormap try *.display()")
        else:
            raise ValueError("ERROR: ds9 cmap = %s, located in %s does not exist" % (ds9_cmap, ds9_cmap_loc))

    def display(self):
        cmap = plt.cm.get_cmap(self.ds9_cmap_name)
        figure = plt.figure(figsize=(14, 10))
        v = np.linspace(0, 1, 4 * cmap.N)

        # Show colormap
        show_cmap = figure.add_axes([0.1, 0.8, 0.8, 0.1])
        im = np.outer(np.ones(50), v)
        show_cmap.imshow(im, cmap=cmap, origin='lower')
        show_cmap.set_xticklabels([])
        show_cmap.set_yticklabels([])
        show_cmap.set_yticks([])
        show_cmap.set_title('RGB colormap of %s' % self.ds9_cmap)

        # Plot RGB profiles
        plot_rgb = figure.add_axes([0.1, 0.1, 0.8, 0.7])
        plot_rgb.plot(v, [cmap(_)[0] for _ in v], color='r')
        plot_rgb.plot(v, [cmap(_)[1] for _ in v], color='g')
        plot_rgb.plot(v, [cmap(_)[2] for _ in v], color='b')
        plot_rgb.set_ylabel('Luminance')
        plot_rgb.set_ylim(-0.005, 1.005)
        plt.show()

    def show(self, image=None, vmin=None, vmax=None, log=True, cmap=None, with_cmap=False, dpi=100, file_reg=None,
             pan_to=None, zoom=None, show_axes=False, figsize=None, savefig=None, fcolor=None, smooth=None):
        if image is None:
            image = 'example.fits'
        if cmap is None:
            cmap = self.ds9_cmap_name
        print("Showing %s" % image)

        hdulist = pyfits.open(image)
        data = hdulist[0].data
        if smooth is not None:
            data=gaussian_filter(data,smooth)
        ny, nx = data.shape
        w = wcs.WCS(hdulist[0].header, relax=False)

        pixim = w.all_world2pix([[10.4099, -9.57457]], 0)
        pixsize = hdulist[0].header['CDELT2'] * 3600.  # pix2asec

        r = np.max(np.array((nx, ny)))
        if zoom != None:
            r = r / zoom

        if pan_to != None:
            if len(pan_to) == 2:
                if zoom == None:
                    pan_to = [pan_to[0], pan_to[1], nx, ny]
                else:
                    pan_to = [pan_to[0], pan_to[1], nx / zoom, ny / zoom]
            elif len(pan_to) == 3:
                pan_to = [pan_to[0], pan_to[1], pan_to[2], pan_to[2]]

            crop_lim = [pan_to[0] - pan_to[2] / 2, pan_to[0] + pan_to[2] / 2, pan_to[1] - pan_to[3] / 2,
                        pan_to[1] + pan_to[3] / 2]
        else:
            if zoom != None:
                pan_to = [nx / 2, ny / 2, nx / zoom, ny / zoom]
                crop_lim = [pan_to[0] - pan_to[2] / 2, pan_to[0] + pan_to[2] / 2, pan_to[1] - pan_to[3] / 2,
                            pan_to[1] + pan_to[3] / 2]

        if figsize is None:
            figsize = (nx / dpi, ny / dpi)

        fig = plt.figure(figsize=figsize, tight_layout=True)
        if with_cmap == True:
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        else:
            ax = fig.add_subplot(111)

        if vmin is None:
            if log == True and np.min(data) <= 0:
                vmin = np.min(data[np.where(data > 0)])
            else:
                vmin = np.min(data)
        if vmax is None:
            vmax = np.max(data)

        if log == True:
            im = ax.imshow(data, cmap=cmap, norm=LogNorm(), vmin=vmin, vmax=vmax, origin='lower')
        else:
            im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')

        if pan_to != None or zoom != None:
            ax.axis(crop_lim)

        if with_cmap == True:
            cax = fig.add_axes([0.95, 0.1, 0.02, 0.8])
            fig.colorbar(im, cax=cax)

        if file_reg is not None:
            for i in range(len(file_reg)):
                draw_reg(ax, file_reg[i], pixsize)

        if fcolor is not None:
            fig.set_facecolor(fcolor)
        else:
            fig.set_facecolor("black")

        if show_axes != True:
            ax.axis('off')

        if savefig is not None:
            fig.savefig(savefig, facecolor=fig.get_facecolor(), edgecolor='none')
        else:
            plt.show()


