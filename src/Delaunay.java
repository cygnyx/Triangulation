// 2d Delaunay Triangulation
// Copyright 1996 (c) Eric C. Olson
// You may use or modify this code for any purpose as long as the
// original author is credited for the original work.

/*
	Algorithm:
	This implementation is derived from a Fortran
	implementation described in a paper for which I don't
	have a reference at the moment.

	The implementation goes as follows:
		create a supertriangle for the unit square.
		scale the points from the bounding box to
			the unit square.
		quick sort the points into scans that sweep
			back and forth across the unit square.
	Now, for each point in the triangulation:
		locate() returns the triangle that the point is in.
		add point by updating edges[] and vertex[] info.
		add new triangles to the testing stack.
	Now, for each testable triangle:
		check point for incircle() enclosure.
		repair edges[] and vertex[]
		add edge triangles to the testing stack.
	Finally, strip off the supertriangle and rescale points.

	The preliminary sort tends to minimize the number of testable
	triangles pushed onto the stack and improves efficiency.  However,
	it isn't necessary for proper triangulation.
*/

package ericco.delaunay;

import java.util.Random;

class Delaunay {
        private static final int MAXSTK = 32;
        private static final int OUTSIDE = 100;
        double  xmin, xmax, ymin, ymax;
        double  x[], y[];
        int     list[];
        double  xrange, yrange, maxrange;
        int     bin[];
        int     location;
        int     stack[];
        int     stktop;
        int     maxstk;
        int     edges[][];
        int     vertex[][];
        int     triangles;

        public double getXMin() { return xmin;}
        public double getYMin() { return ymin;}
        public double getXMax() { return xmax;}
        public double getYMax() { return ymax;}
        public int getNumberPoints() { return list.length - 3;}
        public int getNumberTriangles() { return triangles; }
        public double[] getPoint(int i) {
                double d[] = new double[2];
                d[0] = x[i];
                d[1] = y[i];
                return d;
        }
        public int[] getTriangle(int i) {
                int j[] = new int[3];
                j[0] = vertex[0][i];
                j[1] = vertex[1][i];
                j[2] = vertex[2][i];
                return j;
        }
        public int[] getEdge(int i) {
                int j[] = new int[3];
                j[0] = edges[0][i];
                j[1] = edges[1][i];
                j[2] = edges[2][i];
                return j;
        }

        private int edge(int l, int k) throws Exception {
                if (edges[0][l] == k) return 0;
                if (edges[1][l] == k) return 1;
                if (edges[2][l] == k) return 2;
                throw new Exception("missing edge between "+l+" and "+k);
        }

        private boolean incircle(double x1, double y1, double x2, double y2,
                                 double x3, double y3, double xp, double yp) {
                double x13 = x1 - x3;
                double y13 = y1 - y3;
                double x23 = x2 - x3;
                double y23 = y2 - y3;
                double x1p = x1 - xp;
                double y1p = y1 - yp;
                double x2p = x2 - xp;
                double y2p = y2 - yp;
                double cosa = x13 * x23 + y13 * y23;
                double cosb = x2p * x1p + y1p * y2p;
                if (cosa >= 0.0 && cosb >= 0.0)
                        return false;
                if (cosa <  0.0 && cosb <  0.0)
                        return true;
                double sina = x13 * y23 - x23 * y13;
                double sinb = x2p * y1p - x1p * y2p;
                if (sina * cosb + sinb * cosa < 0.0)
                        return true;
                return false;
        }

        // find triangle containing point
        // vertices are listed in counter clockwise order
        // thus a point inside a triangle must be on the left
        // of all the triangle sides.
        public int locate(double xp, double yp) throws Exception {
                int i = 0, t = location;
                int v1, v2;
//                print(4,"searching for ("+xp+","+yp+")");
                while (i < 3) {
//                        print(5,"checking triangle "+t);
                        v1 = vertex[i][t];
                        v2 = vertex[(i+1)%3][t];
//                        print(5,"current vertex is "+v1);
//                        print(5,"next vertex is "+v2);
                        if ((y[v1]-yp)*(x[v2]-xp) >
                            (x[v1]-xp)*(y[v2]-yp)) {
//                                print(6,"point on right side.");
                                t = edges[i][t];
//                                print(6,"moved to triangle "+t);
                                if (t < 0)
                                        throw new Exception("point off edge");
                                i = 0;
                        } else
                                i++;
                }
//                print(4,"found point in triangle "+t);
                location = t;
                return t;
        }

        private String ptstr(int v) {
                double a = Math.floor(x[v]*100)/100;
                double b = Math.floor(y[v]*100)/100;
                return "(" + a + "," + b + ")";
        }

        private void printtriangle(String pre, int i) {
                String v1, v2, v3;
                int e1, e2, e3;
                v1 = ptstr(vertex[0][i])+vertex[0][i]+" "+edges[0][i];
                v2 = ptstr(vertex[1][i])+vertex[1][i]+" "+edges[1][i];
                v3 = ptstr(vertex[2][i])+vertex[2][i]+" "+edges[2][i];
                System.out.println(pre+i+": "+v1+" "+v2+" "+v3);
        }
/*
        private void dumptriangles () {
                String v1, v2, v3;
                int e1, e2, e3;
                for (int i = 0; i < triangles; i++) {
                        v1 = ptstr(vertex[0][i]);
                        v2 = ptstr(vertex[1][i]);
                        v3 = ptstr(vertex[2][i]);
                        e1 = edges[0][i];
                        e2 = edges[1][i];
                        e3 = edges[2][i];
                        print(5,"tri "+i+": "+v1+" "+e1+" "+v2+" "+e2+" "+
                                              v3+" "+e3);
                }
        }
*/

        private void clearadj(int t) throws Exception {
                int a;
                int i;
                for (i = 0; i < 3; i++) {
                        a = edges[i][t];
                        if (a >= 0)
                                edges[edge(a,t)][a] = -1;
                }
        }

        // strip outside triangles from triangulation
        private void strip() throws Exception {
                int len = list.length - 3;
                int t = 0, s, a;
                while (vertex[0][t] < len &&
                       vertex[1][t] < len &&
                       vertex[2][t] < len) {
//                        printtriangle("keep ",t);
                        t++;
                }

//                printtriangle("drop ",t);
                clearadj(t);

                for (s=t+1; s < triangles; s++) {
                        if (vertex[0][s] >= len ||
                            vertex[1][s] >= len ||
                            vertex[2][s] >= len) {
//                                printtriangle("drop ",s);
                                clearadj(s);
                        } else {
//                                printtriangle("move ",s);
                                for (int i = 0; i < 3; i++) {
                                        a = edges[i][s];
                                        edges[i][t] = a;
                                        vertex[i][t] = vertex[i][s];
                                        if (a >= 0)
                                                edges[edge(a,s)][a] = t;
                                }
                                t++;
                        }
                }
                triangles = t;
        }

        private void triangulate() throws Exception {
                int i, p, l, len = list.length - 3;
                int t;
                double xp, yp;
                int a,b,c;
                int v1, v2, v3;
                int r, erl, era, erb;
//                print(1,"triangulate each point");
                for (i = 0; i < len; i++) {
//                        dumptriangles();
//                        print(2,"Point "+i+" at index "+list[i]);
                        p = list[i];
                        xp = x[p];
                        yp = y[p];
//                        print(3,"Position ("+xp+","+yp+")");
                        t = locate(xp, yp);
//                        print(3,"In triangle "+t);

//                        printtriangle("drop ",t);

                        a = edges[0][t];
                        b = edges[1][t];
                        c = edges[2][t];
                        v1 = vertex[0][t];
                        v2 = vertex[1][t];
                        v3 = vertex[2][t];
                        vertex[0][t] = p;
                        vertex[1][t] = v1;
                        vertex[2][t] = v2;
                        edges[0][t] = triangles + 1;
                        edges[1][t] = a;
                        edges[2][t] = triangles;

//                        printtriangle("add  ",t);

                        vertex[0][triangles] = p;
                        vertex[1][triangles] = v2;
                        vertex[2][triangles] = v3;
                        edges[0][triangles] = t;
                        edges[1][triangles] = b;
                        edges[2][triangles] = triangles + 1;
//                        printtriangle("add  ",triangles);
                        triangles++;

                        vertex[0][triangles] = p;
                        vertex[1][triangles] = v3;
                        vertex[2][triangles] = v1;
                        edges[0][triangles] = triangles - 1;
                        edges[1][triangles] = c;
                        edges[2][triangles] = t;
//                        printtriangle("add  ",triangles);
                        triangles++;

                        if (a >= 0) push(t);
                        if (b >= 0) {
                                edges[edge(b,t)][b] = triangles - 2;
                                push(triangles - 2);
                        }
                        if (c >= 0) {
                                edges[edge(c,t)][c] = triangles - 1;
                                push(triangles - 1);
                        }

                        while (stktop > 0) {
                                l = pop();
                                r = edges[1][l]; // opposite edge
//                                print(6,"comparing "+l+" and "+r);
//                                printtriangle("left ", l);
//                                printtriangle("rght ", r);
                                erl = edge(r,l);
                                era = (erl+1) % 3;
                                erb = (era+1) % 3;
//                                print(6,"lab: "+erl+","+era+","+erb);
                                v1 = vertex[erl][r];
                                v2 = vertex[era][r];
                                v3 = vertex[erb][r];
                                if (incircle(x[v1], y[v1], x[v2], y[v2],
                                             x[v3], y[v3], xp, yp)) {
//                                        print(6,"swapping ...");
                                        a = edges[era][r];
                                        b = edges[erb][r];
                                        c = edges[2][l];

                                        vertex[2][l] = v3;
                                        edges[1][l] = a;
                                        edges[2][l] = r;

                                        vertex[0][r] = p;
                                        vertex[1][r] = v3;
                                        vertex[2][r] = v1;
                                        edges[0][r] = l;
                                        edges[1][r] = b;
                                        edges[2][r] = c;

                                        if (a >= 0) {
                                                edges[edge(a,r)][a] = l;
                                                push(l);
                                        }
                                        if (b >= 0) {
                                                push(r);
                                        }
                                        if (c >= 0) {
                                                edges[edge(c,l)][c] = r;
                                        }
                                }
                        }
                }
                // consistancy check
                if (triangles != 2 * len + 1) {
                        throw new Exception("inconsistent solution");
                }

//                dumptriangles();
//                print(1,"strip external triangles");
                strip();
//                dumptriangles();
        }

        private void push(int item) throws Exception {
                if (stktop == maxstk)
                        throw new Exception("stack overflow");
                stack[stktop++] = item;
//                printtriangle("push ",item);
        }

        private int pop() throws Exception {
                if (stktop == 0)
                        throw new Exception("stack underflow");
//                printtriangle("pop  ",stack[stktop-1]);
                return stack[--stktop];
        }

        public Delaunay(double x[], double y[]) throws Exception {
                int l = Math.min(x.length, y.length);
                if (l < 3) throw new Exception("too few points");
//                print(1, "Triangulate "+l+" points");
//                print(1, "set up "+(l+3)+" points including supertriangle");
                this.x = new double[l+3];
                this.y = new double[l+3];
                this.list = new int[l+3];
//                print(1, "set up "+(2*l+1)+" edges and vertex");
                this.edges = new int[3][2*l+1];
                this.vertex = new int[3][2*l+1];
//                print(1, "set up stack depth to "+l);
                this.stack = new int[l];
                this.maxstk = l;
                this.stktop = 0;
                int v1 = l+0;
                int v2 = l+1;
                int v3 = l+2;
//                print(1, "supertriangle vertex: "+v1+" "+v2+" "+v3);
//                print(1, "vertex 0: set to: "+v1+" "+v2+" "+v3);
                vertex[0][0] = v1;
                vertex[1][0] = v2;
                vertex[2][0] = v3;
//                print(1, "supertriangle edges at 0: set to -1");
                edges[0][0] = -1;
                edges[1][0] = -1;
                edges[2][0] = -1;
                this.x[v1] = -OUTSIDE;
                this.x[v2] =  OUTSIDE;
                this.x[v3] =  0;
                this.y[v1] = -OUTSIDE;
                this.y[v2] = -OUTSIDE;
                this.y[v3] =  OUTSIDE;
                triangles = 1;
//                print(1,"supertriangle vertexs at:");
//                print(2,"("+this.x[v1]+","+this.y[v1]+")");
//                print(2,"("+this.x[v2]+","+this.y[v2]+")");
//                print(2,"("+this.x[v3]+","+this.y[v3]+")");
                location = 0;
//                print(1,"location is 0");
                for (int i = 0; i < l; i++) {
                        this.list[i] = i;
                        xmin = Math.min(x[i],xmin);
                        xmax = Math.max(x[i],xmax);
                        ymin = Math.min(y[i],ymin);
                        ymax = Math.max(y[i],ymax);
                        this.x[i] = x[i];
                        this.y[i] = y[i];
                }
                this.list[l] = l;
                this.list[l+1] = l+1;
                this.list[l+2] = l+2;
//                print(1,"x min, max: "+xmin+","+xmax);
//                print(1,"y min, max: "+ymin+","+ymax);
                xrange = xmax-xmin;
                yrange = ymax-ymin;
                maxrange = Math.max(xrange, yrange);
//                print(1,"range: "+maxrange);
                if (maxrange <= 0.0) throw new Exception("co-located points");
//                print(1,"scale points to unit square");
                for (int i = 0; i < l; i++) {
                        this.x[i] = (this.x[i] - xmin) / maxrange;
                        this.y[i] = (this.y[i] - ymin) / maxrange;
                }

//                dump(2);
//                print(1,"Sort points");
                sort();
//                dump(2);
                triangulate();

//                print(1,"rescale points");
                for (int i = 0; i < l; i++) {
                        this.x[i] = this.x[i] * maxrange + xmin;
                        this.y[i] = this.y[i] * maxrange + ymin;
                }
        }

        private void sort() throws Exception {
                int l = list.length - 3;
                bin = new int[l];
                int partitions = (int)Math.ceil(Math.pow(l, 0.25));
                double xf = partitions / (xrange*1.01/maxrange);
                double yf = partitions / (yrange*1.01/maxrange);
                int i,j,p,k;
                int lstack[] = new int[MAXSTK];
                int rstack[] = new int[MAXSTK];
                int nl, nr;
                int ll = 0;
                int lr = l-1;
                int lm;
                int stktop = 0;
                int guess;
                int t;
                for (k = 0; k < l; k++) {
                        p = list[k];
                        i = (int)(y[p] * yf);
                        j = (int)(x[p] * xf);
                        if (i % 2 == 0)
                                bin[p] = i * partitions + j + 1;
                        else
                                bin[p] = (i+1)*partitions - j;
                }
              for (;;) {
                while (ll < lr) {
                        nl = ll;
                        nr = lr;
                        lm = (ll+lr)/2;
                        guess = bin[list[lm]];
                        for (;;) {
                                while (bin[list[nl]] < guess)
                                        if (++nl == l) break;
                                while (bin[list[nr]] > guess)
                                        if (--nr == 0) break;
                                if (nl < nr-1) {
                                        t = list[nl];
                                        list[nl] = list[nr];
                                        list[nr] = t;
                                        nl++;
                                        nr--;
                                } else
                                        break;
                        }
                        if (nl <= nr) {
                                if (nl < nr) {
                                        t = list[nl];
                                        list[nl] = list[nr];
                                        list[nr] = t;
                                }
                                nl++;
                                nr--;
                        }
                        if (nr < lm) {
                                lstack[stktop] = nl;
                                rstack[stktop] = lr;
                                lr = nr;
                        } else {
                                lstack[stktop] = ll;
                                rstack[stktop] = nr;
                                ll = nl;
                        }
                        if (++stktop == MAXSTK)
                                throw new Exception("stack overflow");
                }
                if (stktop == 0) break;
                stktop--;
                ll=lstack[stktop];
                lr=rstack[stktop];
              }
        }

        private double[] center(int v1, int v2, int v3) {
                double x1 = x[v1];
                double y1 = y[v1];
                double x2 = x[v2];
                double y2 = y[v2];
                double x3 = x[v3];
                double y3 = y[v3];
                double m[] = new double[2];
                double x12 = (x1 + x2) / 2;
                double y12 = (y1 + y2) / 2;
                double x23 = (x2 + x3) / 2;
                double y23 = (y2 + y3) / 2;
                double dx12 = (x1 - x2);
                double dy12 = (y1 - y2);
                double dx23 = (x2 - x3);
                double dy23 = (y2 - y3);
                double d = dy12 * dx23 - dy23 * dx12;
                double x = dy12 * dy23 * (y23 - y12);
                x = x + dy12 * dx23 * x23;
                x = x - dy23 * dx12 * x12;
                try {
                        m[0] = x / d;
                } catch (Exception e) {
                        m[0] = (x >= 0) ? 1.6e38 : -1.6e38;
                }
                try {
                        m[1] = y23 - dx23 * (m[0] - x23) / dy23;
                } catch (Exception e) {
                        try {
                        m[1] = y12 - dx12 * (m[0] - x12) / dy12;
                        } catch (Exception e1) {
                                m[1] = (y1 > y3) ? 1.6e38 : -1.6e38;
                        }
                }
                return m;
        }

        private double distance(double x1, double y1, double x2, double y2) {
                double dx = x1 - x2;
                double dy = y1 - y2;
                return Math.sqrt(dx*dx+dy*dy);
        }

        private void bruteforcecheck() throws Exception {
                int i, j, len = list.length - 3;
                int v1, v2, v3;
                double c[], cx, cy, r;
                for (i = 0; i < triangles; i++) {
                        v1 = vertex[0][i];
                        v2 = vertex[1][i];
                        v3 = vertex[2][i];
                        c = center(v1,v2,v3);
                        cx = c[0];
                        cy = c[1];
                        r = distance(cx, cy, x[v1], y[v1]);
                        for (j = 0; j < len; j++) {
                                if (j == v1 || j == v2 || j == v3)
                                        continue;
                                if (r > 1.01 * distance(cx, cy, x[j], y[j]))
                                        throw new Exception(
 "not Delaunay: "+ptstr(j)+" in triangle: "+
  ptstr(v1)+","+ptstr(v2)+","+ptstr(v3));
                        }
                }
        }

/*
        private void dump(int indent) {
                int p, l = list.length - 3;
                if (bin != null)
                for (int i = 0; i < l; i++) {
                        p = list[i];
                        print(indent,
                                "Point "+p+" "+bin[p]+" ("+x[p]+","+y[p]+")");
                }
                else
                for (int i = 0; i < l; i++) {
                        p = list[i];
                        print(indent,
                                "Point "+p+" ("+x[p]+","+y[p]+")");
                }
        }
*/

        private static void print(int i, String str) {
                while (i-- != 0) System.out.print(" ");
                System.out.println(str);
        }

        final static int DEFAULT_POINTS = 20;

	public static void main (String args[]) {
                print(0,"Delaunay Triangulation Test");

                Random  rnd = new Random(10000);
                int     points = DEFAULT_POINTS;
                double  x[] = new double[points];
                double  y[] = new double[points];

                print(0,"Generate random points");
                for (int i = 0; i < points; i++) {
                        x[i] = rnd.nextDouble();
                        y[i] = rnd.nextDouble();
                        print(1,"Point "+i+" ("+x[i]+","+y[i]+")");
                }

                try {
//                print(0,"Create Triangulation from random points");
                Delaunay d = new Delaunay(x,y);
//                print(0,"Delaunay Triangulation completed");
                d.bruteforcecheck();
                } catch (Exception e) {
                        print(0,e.toString());
                }
		System.exit(0);
	}
}
