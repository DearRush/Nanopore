# -*- coding: utf-8 -*-
import numpy as np
import matplotlib as plt
from NDFs import sNDF
ABC = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T']

ax = []

x = np.linspace(-5, 5, 101)


fig_size = ( 3.25*len( Bx ), 2 )
fig = plt.figure( figsize=fig_size )
gs = plt.GridSpec( 1, len( Bx ), top=0.95, wspace=1, hspace=0.75)

for c, B in enumerate( Bx ):
    ax.append( fig.add_subplot( gs[0, c] ) )
    
    HWHM = S * np.sqrt( ( 2 * np.log(2)**(1/B) ) )
    low_HWHM = ( ( x0 - min(x) ) - HWHM ) / ( max(x)-min(x) )
    
    ax[-1].axhline( y=-5, xmin=low_HWHM, xmax=1-low_HWHM, linestyle='-', color='k', linewidth=2 )
    ax[-1].axvline( x=HWHM, ymin=(20/130)-0.02, ymax=(20/130)+0.02, color='k', linewidth=2 )
    ax[-1].axvline( x=-HWHM, ymin=(20/130)-0.02, ymax=(20/130)+0.02, color='k', linewidth=2 )
    ax[-1].axvline( x=HWHM, ymin=(20/130)+0.02, ymax=((((a0*100)/2)+25)/130), color='k', linewidth=2, alpha=0.5, linestyle=':' )
    ax[-1].axvline( x=-HWHM, ymin=(20/130)+0.02, ymax=((((a0*100)/2)+25)/130), color='k', linewidth=2, alpha=0.5, linestyle=':' )
    
    ax[-1].plot( x, 100*(sNDF( x, a0, x0, S, B, Offset )+((np.random.rand(len(x))-0.5)*0.1)), 'k', linewidth=2, alpha=0.5 )
    ax[-1].plot( x, 100*(sNDF( x, a0, x0, S, B, Offset )), 'r', linewidth=1.5 )
    ax[-1].axvline( x=0, ymin=25/130, linestyle=':', color='k', linewidth=2, alpha=0.5 )
    
    ax[-1].axhline( a0*100, xmin=0.5, color='r', linestyle=':', linewidth=2 )
    ax[-1].axhline( Offset*100, color='r', linestyle=':', linewidth=2 )
    
    ax[-1].axvline( x=4.5, ymin= 25/130, ymax=(a0*100+25)/130, color='k', linewidth=2 )
    ax[-1].axhline( y=a0*100, xmin=0.95-0.02, xmax=0.95+0.02, color='k', linewidth=2 )
    ax[-1].axhline( y=Offset, xmin=0.95-0.02, xmax=0.95+0.02, color='k', linewidth=2 )
    
    ax[-1].set_ylabel('$I_{excluded}$ (%)', fontsize=fontsize)
    ax[-1].set_ylim([-25, 105])
    ax[-1].set_xlim([-5, 5])
    ax[-1].text(5.3, 80, '$I_{B}$', fontsize=12, weight='bold', va='center')
    ax[-1].text(0, -21, '$FWHM$', horizontalalignment='center', fontsize=12, weight='bold', ha='center')
    ax[-1].text(5.3, 0, '$I_{o}$', fontsize=12, weight='bold', va='center')
    ax[-1].text(-4.5, 90, r'$\beta$ = ' + str( B ), fontsize=12, weight='bold', va='center')
    ax[-1].text(5.3, 40, r'$\Delta I_{B}$', fontsize=12, weight='bold', va='center')
    ax[-1].text(0, 110, r'$\mu$', fontsize=12, weight='bold', ha='center')
    
    ax[-1].set_xticks([])
    ax[-1].text(-9.2, 120, ABC[c], fontsize=24, weight='bold')

