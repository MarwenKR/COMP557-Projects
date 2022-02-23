/*
 * Marwen Kraiem 260955202
 */

package comp557.a4;

import javax.vecmath.Color4f;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

public class Light {
	
	/** Light name */
    public String name = "";
    
    /** Light colour, default is white */
    public Color4f color = new Color4f(1,1,1,1);
    
    /** Light position, default is the origin */ //will be used as a corner if the type of light is area
    public Point3d from = new Point3d(0,0,0);
    
    /** Light intensity, I, combined with colour is used in shading */
    public double power = 1.0;
    
    /** Type of light, default is a point light */
    public String type = "point";

    
    /** light samples */
    public int samples = 1;
    public Vector3d u = new Vector3d(0,0,1);
    public Vector3d v = new Vector3d(1,0,0);
    
    /**
     * Default constructor 
     */
    public Light() {
    	// do nothing
    }
}
