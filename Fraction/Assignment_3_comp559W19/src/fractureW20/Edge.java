package fractureW20;

import javax.vecmath.Point2d;

/**
 * Edge class to help with collision detection, and drawing the current boundary.
 * We store the triangle opposite the edge to help in collision detection, but
 * note that the list of all adjacent triangles to a particle is in the particle!
 * @author kry
 */
public class Edge {

    Particle p1 = null;
    
    Particle p2 = null;
    
    /** triangle adjacent to this edge */
    FEMTriangle triangle;
        
    /**
     * Creates an Edge connecting two particles.
     * The rest length should be set
     * @param p1
     * @param p2
     */
    public Edge( Particle p1, Particle p2, FEMTriangle triangle ) {
        this.p1 = p1;
        this.p2 = p2;
        this.triangle = triangle;        
    }
    
    /** 
     * Checks to see if this edge goes through the same points, disregarding order
     * (though in some scenarios we may care a lot about the order 
     */
    @Override
    public boolean equals(Object obj) {
    	Edge e = (Edge) obj;
    	return (p1 == e.p1 && p2 == e.p2) || (p1 == e.p2 && p2 == e.p1);
    }
    
    @Override
    public int hashCode() {
    	final Point2d p = new Point2d();
    	p.add(p1.p,p2.p);
    	return p.toString().hashCode(); // gross!  but works.  	
    }
    
    
}
