# Code to bypass ds9 and use python to display images


 first load the ds9 colormap, e.g.:

 cm=ds9_colormap(CIAO='/home/vittorio/Software/ciao-4.11',ds9_cmap='sls')

 then just display an image using e.g.

 cm.show(image='example.fits')

 or with preferred options:

 cm.show(image='example.fits',file_reg='test.reg',pan_to=[420,430],zoom=1,show_axes=False,with_cmap=True, vmin=7e-7, vmax=1e-3, log=True, cmap='ds9sls')

 Notes:

 an_to   : center the image around this point, if a third/forth value is given, then also crop image with a rectangle around he point using given radius
 zoom     : self-explanatory
 show_axes: to display the axes on the image 
 cmap     : can use different cmap
 with_cmap: to display colorbar

