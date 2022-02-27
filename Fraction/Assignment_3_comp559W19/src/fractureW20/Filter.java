package fractureW20;

import no.uib.cipr.matrix.Vector;

public interface Filter {

    /**
     * Removes disallowed parts of v by projection
     * @param v
     */
    public void filter( Vector v );
    
}
