package fractureW20;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector2d;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

/**
 * Spring between a FEMTriangle and a point driven by the mouse position
 * @author kry
 */
public class MouseSpring {

    FEMTriangle triangle;
    
    Point3d barycentricCoords = new Point3d();
    
    Point2d mousePos = new Point2d();
        
    double k = 100;

    double b = 1;
    
    /** 
     * Nothing to initialize as this will be set later 
     */ 
    public MouseSpring() {
    	// do nothing
    }

    public void setMousePosition( Point2d p ) {
    	mousePos.set( p );
    }
    
    public void setFEMTriangle( FEMTriangle t, Point3d coords ) {
    	this.triangle = t;
    	this.barycentricCoords.set( coords );
    }
    
    /**
     * Applies explicit spring forces to the triangle's nodes.
     */
    public void apply() {
    	if ( triangle == null ) return;
    	
    	final Point2d p = new Point2d();
    	final Vector2d tmp = new Vector2d();
    	tmp.scale( barycentricCoords.x, triangle.A.p );
    	p.add(tmp);
    	tmp.scale( barycentricCoords.y, triangle.B.p );
    	p.add(tmp);
    	tmp.scale( barycentricCoords.z, triangle.C.p );
    	p.add(tmp);
    	
        final Vector2d force = new Vector2d();
    	force.sub( mousePos, p );
        force.scale( k );
        
        // again use barycentric coordinates to distribute the force
        tmp.scale( barycentricCoords.x, force );
        triangle.A.addForce(tmp);
        tmp.scale( barycentricCoords.y, force );
        triangle.B.addForce(tmp);
        tmp.scale( barycentricCoords.z, force );
        triangle.C.addForce(tmp);
        
    }
    
    public void display( GLAutoDrawable drawable ) {
    	if ( triangle == null ) return;
    	GL2 gl = drawable.getGL().getGL2();
    	final Point2d p = new Point2d();
    	final Vector2d tmp = new Vector2d();
    	tmp.scale( barycentricCoords.x, triangle.A.p );
    	p.add(tmp);
    	tmp.scale( barycentricCoords.y, triangle.B.p );
    	p.add(tmp);
    	tmp.scale( barycentricCoords.z, triangle.C.p );
    	p.add(tmp);
	    gl.glColor4d( 1,0,0, 1 );
		gl.glBegin( GL.GL_LINES );
	    gl.glVertex2d( mousePos.x, mousePos.y );
	    gl.glVertex2d( p.x, p.y );
	    gl.glEnd();
    }
    
}
