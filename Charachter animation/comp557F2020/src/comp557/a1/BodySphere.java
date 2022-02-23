package comp557.a1;

import javax.vecmath.Tuple3d;

import com.jogamp.opengl.GL4;
import com.jogamp.opengl.GLAutoDrawable;

import comp557.a1.geom.Sphere;

//**** Marwen Kraiem 260955202 ****//

public class BodySphere extends GraphNode {


	//dimensions
	double cx,cy,cz;
	
	//angles
	double rx,ry,rz;
	
	//scale
	double sx,sy,sz;
	
	//color
	float r, g, b;
		
	public BodySphere( String name, double cx, double cy, double cz, double sx, double sy, double sz, double rx, double ry, double rz, double r, double g, double b ) {
		super(name);
		
		// center
		this.cx = cx;
		this.cy = cy;
		this.cz = cz;
		
		// rotation
		this.rx = rx;
		this.ry = ry;
		this.rz = rz;
		
		// scale
		this.sx = sx;
		this.sy = sy;
		this.sz = sz;
		
		// color
		this.r = (float) r;
		this.g = (float) g;
		this.b = (float) b;


	}

	public void setCentre(Tuple3d t  ) {
		cx = t.x;
		cy = t.y;
		cz = t.z;

	}
	
	public void setScale( Tuple3d t) {
		sx = t.x;
		sy = t.y;
		sz = t.z;
	}
	
	public void setRotation( Tuple3d t) {
		rx = t.x;
		ry = t.y;
		rz = t.z;
	}
	public void setColor(Tuple3d t  ) {
		r = (float) t.x;
		g = (float) t.y;
		b = (float) t.z;
	}
	
	@Override
	public void display( GLAutoDrawable drawable, BasicPipeline pipeline ) {
		pipeline.push();
		GL4 gl = drawable.getGL().getGL4();		
		gl.glUniform3f( pipeline.kdID, r, g, b );
		gl.glUniform3f( pipeline.kaID, r, g, b );
		gl.glUniform3f( pipeline.ksID,0.5f, 0.5f,0.5f );
        gl.glUniform1f( pipeline.ShininessID, 32.00f );

		pipeline.translate(cx,cy,cz);
		pipeline.scale(sx,sy,sz);
		pipeline.rotate(rx,1,0,0);
		pipeline.rotate(ry,0,1,0);
		pipeline.rotate(rz,0,0,1);

		Sphere.draw(drawable, pipeline);
		pipeline.pop();
	}
	
}
