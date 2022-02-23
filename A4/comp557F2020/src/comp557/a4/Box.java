/*
 * Marwen Kraiem 260955202
 */
package comp557.a4;

import javax.vecmath.Point3d;

/**
 * A simple box class. A box is defined by it's lower (@see min) and upper (@see max) corner. 
 */
public class Box extends Intersectable {

	public Point3d max;
	public Point3d min;
	
    /**
     * Default constructor. Creates a 2x2x2 box centered at (0,0,0)
     */
    public Box() {
    	super();
    	this.max = new Point3d( 1, 1, 1 );
    	this.min = new Point3d( -1, -1, -1 );
    }	

	@Override
	public void intersect(Ray ray, IntersectResult result) {
		// TODO: Objective 6: intersection of Ray with axis aligned box
		double tmin , tmax, txmin, txmax, tymin, tymax, tzmin, tzmax;
			txmin = Math.abs(ray.viewDirection.x)>1e-9? Math.min((min.x-ray.eyePoint.x)/ray.viewDirection.x,(max.x-ray.eyePoint.x)/ray.viewDirection.x): Math.min((min.x-ray.eyePoint.x),(max.x-ray.eyePoint.x))*Math.signum(ray.viewDirection.x)*Double.POSITIVE_INFINITY;
			txmax = Math.abs(ray.viewDirection.x)>1e-9? Math.max((min.x-ray.eyePoint.x)/ray.viewDirection.x,(max.x-ray.eyePoint.x)/ray.viewDirection.x): Math.max((min.x-ray.eyePoint.x),(max.x-ray.eyePoint.x))*Math.signum(ray.viewDirection.x)*Double.POSITIVE_INFINITY;
	
			tymin = Math.abs(ray.viewDirection.y)>1e-9? Math.min((min.y-ray.eyePoint.y)/ray.viewDirection.y,(max.y-ray.eyePoint.y)/ray.viewDirection.y):Math.min((min.y-ray.eyePoint.y),(max.y-ray.eyePoint.y))*Math.signum(ray.viewDirection.y)*Double.POSITIVE_INFINITY;	
			tymax = Math.abs(ray.viewDirection.y)>1e-9? Math.max((min.y-ray.eyePoint.y)/ray.viewDirection.y,(max.y-ray.eyePoint.y)/ray.viewDirection.y): Math.max((min.y-ray.eyePoint.y),(max.y-ray.eyePoint.y))*Math.signum(ray.viewDirection.y)*Double.POSITIVE_INFINITY;		
	
			tzmin = Math.abs(ray.viewDirection.z)>1e-9? Math.min((min.z-ray.eyePoint.z)/ray.viewDirection.z,(max.z-ray.eyePoint.z)/ray.viewDirection.z): Math.min((min.z-ray.eyePoint.z),(max.z-ray.eyePoint.z))*Math.signum(ray.viewDirection.z)*Double.POSITIVE_INFINITY;
			tzmax = Math.abs(ray.viewDirection.z)>1e-9? Math.max((min.z-ray.eyePoint.z)/ray.viewDirection.z,(max.z-ray.eyePoint.z)/ray.viewDirection.z): Math.max((min.z-ray.eyePoint.z),(max.z-ray.eyePoint.z))*Math.signum(ray.viewDirection.z)*Double.POSITIVE_INFINITY;

		
		tmin = Math.max(txmin, Math.max(tymin, tzmin));
		tmax = Math.min(txmax, Math.min(tymax, tzmax));

		if(tmin<tmax && tmin>1e-9 && tmin < result.t) {
			
			// intersection point
			result.t = tmin;		
		    result.p.scaleAdd(tmin, ray.viewDirection, ray.eyePoint);
		   
		    // material 
		    result.material = material;
		    
		    // normal vector 
		    if(Math.abs(result.p.x-min.x)<1e-9)
				result.n.set(-1,0,0);
			else if(Math.abs(result.p.x-max.x)<1e-9)
				result.n.set(1,0,0);
			else if(Math.abs(result.p.y-min.y)<1e-9)
				result.n.set(0,-1,0);
			else if(Math.abs(result.p.y-max.y)<1e-9)
				result.n.set(0,1,0);
			else if(Math.abs(result.p.z-min.z)<1e-9)
				result.n.set(0,0,-1);
			else if(Math.abs(result.p.z-max.z)<1e-9)
				result.n.set(0,0,1);
			else result.n.set(1,1,1);
			
		      
		}

	
		
	}	

}
