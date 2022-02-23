/*
 * Marwen Kraiem 260955202
 */

package comp557.a4;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;


public class Quadric extends Intersectable {
    
	/**
	 * Quadric elements.
	 */
	public Matrix4d Q = new Matrix4d();
	public Matrix3d A = new Matrix3d();
	public Vector3d B = new Vector3d();
	public double C;
	
	/**
	 * The second material, e.g., for front and back?
	 */
	Material material2 = null;
   
	public Quadric() {
    	
	}
	
	@Override
	public void intersect(Ray ray, IntersectResult result) {
		
		Vector3d Ad = new Vector3d();
		A.transform(ray.viewDirection, Ad);
		Vector3d Ap = new Vector3d();
		A.transform(ray.eyePoint, Ap);

		Vector3d p = new Vector3d(ray.eyePoint.x,ray.eyePoint.y,ray.eyePoint.z);
		double a= ray.viewDirection.dot(Ad);
		double b = ray.viewDirection.dot(Ap)- B.dot(p);
		double c = C-2*B.dot(p)+p.dot(Ap);
		
		double delta = b*b-a*c;
		if (delta >= 0) {	
    		double t = (-b - Math.sqrt(delta))/a;

	        if (  t > 1e-9 && t < result.t) {
	        	 result.t = t;
	             result.material = material;
	             result.p.scaleAdd(t, ray.viewDirection, ray.eyePoint);
	             A.transform(result.p,result.n);	  
	             result.n.sub(B);	          
	             result.n.normalize();
	             
	        	 
	        }
		}

	}
	
}
