      subroutine find_wall_2D(on_wall,dt,iface,r2i,r2ip1,
     &     tan_th2j,tan_th2jp1,ux,uy,uz,x1,y1,z1)
c
      use closest_wall
  
c   This routine calculates the distance to all faces of a 2D grid cell
c   (in a spherical -- r,theta) from a point (x1,y1,z1) along the
c   directional cosine vector (ux,uy,yz).
c     NOTE: ***** it is necessary to flag tan_th2 array using a -1.
c               to indicate theta=90.
c        Program ASSUMES that theta=90 can occur in ONLY ONE of the
c        tan_th2 values.
c
c  Output values:
c     dt = distance from (x1,y1,z1) along (ux,uy,uz) to closest face(s),
c          along for possibility of passing through cell vertex (i.e.,
c     ioffset = integer offsets for new cell position.  
c          intersecting two constant surfaces).
c
c  Input values:
c     r2i,r2ip1 = squared radii of two constant surfaces r(i),r(i+1)
c     tan_th2j,jp1 = (tan(theta))**2 for two constant surfaces
c        theta(j),theta(j+1).  see NOTE above.
c     ux,uy,uz = directional cosines along x,y,z directions
c     x1,y1,z1 = point from which distance is calculated
c
c m.j. wolff/b.a. whitney, 2000/03/16 (created from 3D version)
c history:
c


      implicit none
c ... 
      real*8 dt,r2i,r2ip1,tan_th2j,tan_th2jp1,ux,uy,uz,x1,y1,z1
      logical(1) :: on_wall
c ..
      integer iface
      
      call reset_t()

c find shortest non-negative distance to two constant radius surfaces
c where negative ifound = .false. indicates no intersection    
      call radface(ux,uy,uz,x1,y1,z1,r2i,r2ip1)

c find shortest positive distance to two constant theta surfaces
c where negative ifound = .false. indicates no intersection    
      call thetaface(ux,uy,uz,x1,y1,z1,tan_th2j,tan_th2jp1)

      call find_next_wall(on_wall,dt,iface)

      return
      end
