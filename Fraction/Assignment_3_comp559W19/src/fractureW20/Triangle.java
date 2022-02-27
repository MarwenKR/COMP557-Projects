package fractureW20;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;

import javax.vecmath.Color3f;
import javax.vecmath.Matrix3d;
import javax.vecmath.Vector3d;

public class Triangle {

	public static void main(String args[]) {
		//runTriangle("triangle --help");
		String fnameroot = "dataFracture/square";
        genTriangles( fnameroot );
		runTriangle( "triangle -pa500qD " + fnameroot + ".poly" );
    }
	
	static void runTriangle( String command ) {
		String s = null;
        try {    
        	Process p = Runtime.getRuntime().exec( command );
            BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
            BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            // read the output from the command
            System.out.println("Here is the standard output of the command:\n");
            while ((s = stdInput.readLine()) != null) {
                System.out.println(s);
            }
            // read any errors from the attempted command
            System.out.println("Here is the standard error of the command (if any):\n");
            while ((s = stdError.readLine()) != null) {
                System.err.println(s);
            }
        } catch (IOException e) {
            System.out.println("exception happened - here's what I know: ");
            e.printStackTrace();
        }
	}
	
	static void genTriangles( String fnameroot) {
		try {
			FileOutputStream fos = new FileOutputStream( new File(fnameroot + ".poly") );
			PrintStream ps = new PrintStream( fos );
			
			int n = 4;
			double r = 70;
			double[] x = new double[n];
			double[] y = new double[n];
			for ( int i = 0; i < n; i++ ) {
				double theta = (i / (double) n) * Math.PI * 2 + Math.PI / 4.0;
				x[i] = Math.cos(theta) * r;
				y[i] = Math.sin(theta) * r;
			}
			// https://www.cs.cmu.edu/~quake/triangle.poly.html
			// We'll use dim 2, attributes 2, markers 0, and set
			// the attributes to be suitable texture coordinates 
			ps.println( n + " 2 0 0" ); // vertices
			for ( int i = 0; i < n; i++ ) {
				ps.println( (i+1) + " " + x[i] +" " + y[i] ); // + " " + x[i] +" " + y[i] );
			}
			ps.println( n + " 0" ); // segments
			for ( int i = 0; i < n; i++ ) {
				ps.println( (i+1) + "  " + (i+1) + " " + (((i+1)%n)+1) );
			}
			ps.println("0"); // holes
			ps.close();
			fos.close();
		} catch ( Exception e ) {
			e.printStackTrace();
		}
	}
	
	static void loadTriangles( FEMSystem system, String fnameroot, Matrix3d M ) {
		try {
			// Node file
			// First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
			// Remaining lines: <vertex #> <x> <y> [attributes] [boundary marker]
			FileInputStream fis = new FileInputStream( new File(fnameroot + ".node") );
			InputStreamReader isr = new InputStreamReader( fis );
			BufferedReader reader = new BufferedReader( isr );
			Vector3d v = new Vector3d();
			Particle[] particles = null; // particles, for indices in loading the elements
			while (true) {
				String line = reader.readLine();
				if ( line == null ) break;
				if ( line.charAt(0) == '#') continue;
				line = line.trim();
				String[] fields = line.split("\\s+");
				int n = Integer.parseInt(fields[0]);  // number of nodes
				/*int na = */ Integer.parseInt(fields[2]); // number of attributes, not needed		
				particles = new Particle[n]; // to keep track of indices for triangles
				int count = 0;
				while (true) {
					line = reader.readLine();
					if (line == null) break;
					if ( line.charAt(0) == '#') continue;
					line = line.trim();
					String[] data = line.split("\\s+");
					v.set( Double.parseDouble(data[1]), Double.parseDouble(data[2]), 1 );
					M.transform(v);
					Particle p = system.createParticle( v.x, v.y, 0 , 0 );
					p.mass = 0; // we'll set this later based on triangle area
					particles[count++] = p;
					if ( count >= n ) break;
				}
			}
			reader.close();
			isr.close();
			fis.close();
			// Element file
			// First line: <# of triangles> <nodes per triangle> <# of attributes>
			// Remaining lines: <triangle #> <node> <node> <node> ... [attributes]
			fis = new FileInputStream( new File(fnameroot + ".ele") );
			isr = new InputStreamReader( fis );
			reader = new BufferedReader( isr );
			double density = 1e-2; // no real units here sadly
			Color3f colour = new Color3f(0,1,0);
			while (true) {
				String line = reader.readLine();				
				if ( line == null ) break;
				if ( line.charAt(0) == '#') continue;
				line = line.trim();
				String[] fields = line.split("\\s+");
				int n = Integer.parseInt(fields[0]);  // number of triangles
				int count = 0;
				while (true) {
					line = reader.readLine();
					if (line == null) break;
					if ( line.charAt(0) == '#') continue;
					line = line.trim();
					String[] data = line.split("\\s+");
					int i1 = Integer.parseInt(data[1])-1;
					int i2 = Integer.parseInt(data[2])-1;
					int i3 = Integer.parseInt(data[3])-1;
					FEMTriangle tri = new FEMTriangle( particles[i1], particles[i2], particles[i3], colour );
					tri.density = density;
					particles[i1].mass += 1.0/3.0 * tri.area * density;
					particles[i2].mass += 1.0/3.0 * tri.area * density;
					particles[i3].mass += 1.0/3.0 * tri.area * density;
					system.femSprings.add( tri );
					count++;
					if ( count >= n ) break;
				}
			}
			reader.close();
			isr.close();
			fis.close();

			// we'll find the edges based on the elements that were loaded, much like the way images are loaded.
			
		} catch ( Exception e ) {
			e.printStackTrace();
		}	
	}
	
}
