// Triangulation Viewer
// Copyright 1996 (c) Eric C. Olson
// You may use or modify this code for any purpose as long as the
// original author is credited for the original work.

/*
	My only excuse is that this was one of my first attempts
	to write Java code.  Many improvements are needed, including:
		renaming of gui widgets and classes.
		try to follow a naming convention.
		cleanup of duplicate code.
		modular organization of code.
		documentation.
		string dictionary for translations.
		fixup callbacks.
		use awt Point.
		subclass Polygon to triangles.
		do something with longLine.
		fix up color usage.
		perhaps remove XOR mode ...
		remove complex cases in paint() 
*/

package ericco.delaunay;

import java.applet.Applet;

import java.awt.*;

import java.util.Random;


public class TriangulationViewer extends Applet {
public String getAppletInfo () {
        return "TriangulationViewer v0.1 by Eric C. Olson";
}
public String[][] getParameterInfo(){
        String p[][] = {
                {"data",         "string",       "data set to plot"},
                {"label",        "string",       "plot label"},
                {"points",       "int",          "number of random points"},
                {"display",      "string",       "display format"},
                {"foreground",   "string",       "foreground format"},
                {"background",   "string",       "background color"},
        };
        return p;
}

Choice  display;
TriangulationCanvas tc;
Random  rnd = null;
private String displaytype = null;
private boolean pointsset = false;
protected int numpoints = 20;
private String datapoints = null;
private String label = null;
private String foreground = null;

public void init() {
        showStatus("Loading TriangulationViewer ...");
        displaytype = getParameter("display");
        datapoints = getParameter("data");
        label = getParameter("label");
        foreground = getParameter("foreground");

        String background = getParameter("background");
        Color bgcolor = null;
	if (background != null) bgcolor = Color.getColor(background);
        if (bgcolor != null) setBackground(bgcolor);

        String p = getParameter("points");
        if (p != null) {
                pointsset = true;
                numpoints = Integer.valueOf(p).intValue();
        }
        try {
                Delaunay t;
                if (datapoints == null)
                        t = RandomPointSet(numpoints);
                else
                        t = ParsePointSet(datapoints);
                tc = new TriangulationCanvas(this,t);
                if (displaytype != null)
                        setdisplay(displaytype);
                if (foreground != null)
                        setforeground(foreground);
                setLayout(new BorderLayout());
                Panel north = northpanel();
                if (north != null)
                        add("North", north);
                add("Center", tc);
                if (label != null)
                        add("South", new Label(label));
        } catch (Exception e) {
                showStatus("unloadable: "+e.toString());
        }
}

private Delaunay ParsePointSet(String points) throws Exception {
        int     maxlen = points.length()/4 + 1;
        int     i = 0;
        int     lastindex;
        int     index, q;
        double  x[], mx[] = new double[maxlen];
        double  y[], my[] = new double[maxlen];

        for (lastindex = 0; true; lastindex = index + 1) {
                index = points.indexOf(',', lastindex);
                if (index == -1)
                        break;
                mx[i] = Double.valueOf(points.substring(lastindex, index)
                                ).doubleValue();
                lastindex = index + 1;
                q = index = points.indexOf(',', lastindex);
                if (q == -1) q = points.length();
                my[i] = Double.valueOf(points.substring(lastindex, q)
                                ).doubleValue();
                i++;
                if (index == -1)
                        break;
        }
        x = new double[i];
        y = new double[i];
        for (int j = 0 ; j < i; j++) {
                x[j] = mx[j];
                y[j] = my[j];
        }
        return new Delaunay(x, y);
}
private Delaunay RandomPointSet(int points) throws Exception {
        if (rnd == null) rnd = new Random(10000);
        double  x[] = new double[points];
        double  y[] = new double[points];

        for (int i = 0; i < points; i++) {
                x[i] = rnd.nextDouble();
                y[i] = rnd.nextDouble();
        }
        return new Delaunay(x, y);
}

private Panel northpanel() {
        Panel p = null;
        if (datapoints == null) {
                if (p == null) p = new Panel();
                p.add(new NewButton(this));
        }

        if (displaytype == null) {
                if (p == null) p = new Panel();
                p.add(new Label("Display:"));
                p.add(new DisplayChoice(this));
        }

        if (foreground == null) {
                if (p == null) p = new Panel();
                p.add(new Label("Highlight:"));
                p.add(new HighlightChoice(this));
        }
        if (datapoints == null && pointsset == false) {
                if (p == null) p = new Panel();
                p.add(new Label("Points:"));
                p.add(new PointChoice(this));
        }

        return p;
}


protected void setdisplay(String s) {
        int d = TriangulationCanvas.POINTS;
        if (s.equals("Points")) {
                d = TriangulationCanvas.POINTS;
        } else if (s.equals("Triangles")) {
                d = TriangulationCanvas.TRIANGLES;
        } else if (s.equals("Circles")) {
                d = TriangulationCanvas.CIRCLES;
        } else if (s.equals("Centers")) {
                d = TriangulationCanvas.CENTERS;
        } else if (s.equals("Border")) {
                d = TriangulationCanvas.BORDER;
        } else if (s.equals("Tessellation")) {
                d = TriangulationCanvas.TESSELLATION;
        }
        tc.setDisplay(d);
}

protected void setforeground(String s) {
        int d = TriangulationCanvas.NODRAW;
        if (s.equals("fillOval")) {
                d = TriangulationCanvas.FILLOVAL;
        } else if (s.equals("drawOval")) {
                d = TriangulationCanvas.DRAWOVAL;
        } else if (s.equals("drawTriangle")) {
                d = TriangulationCanvas.DRAWTRIANGLE;
        } else if (s.equals("fillTriangle")) {
                d = TriangulationCanvas.FILLTRIANGLE;
        }
        tc.setforeground(d);
}

protected void newrandom() {
        try {
                Delaunay t = RandomPointSet(numpoints);
                tc.setTriangulation(t);
        } catch (Exception e) {
                showStatus("couldn't restart");
        }
}
}

class NewButton extends Button {

TriangulationViewer tv;

NewButton(TriangulationViewer tv) {
        super("New Data Set");
        this.tv = tv;
}

public boolean action (Event evt, Object arg) {
        tv.newrandom();
        return true;
}
}

class DisplayChoice extends Choice {

TriangulationViewer tv;

DisplayChoice(TriangulationViewer tv) {
        super();
        this.tv = tv;
        addItem("Points");
        addItem("Triangles");
        addItem("Circles");
        addItem("Centers");
        addItem("Border");
        addItem("Tessellation");
        select(1);
        tv.setdisplay(getSelectedItem());
}

public boolean action (Event evt, Object arg) {
        tv.setdisplay(getSelectedItem());
        return true;
}
}

class HighlightChoice extends Choice {

TriangulationViewer tv;

HighlightChoice(TriangulationViewer tv) {
        super();
        this.tv = tv;
        addItem("none");
        addItem("drawOval");
        addItem("fillOval");
        addItem("drawTriangle");
        addItem("fillTriangle");
        select(1);
        tv.setforeground(getSelectedItem());
}

public boolean action (Event evt, Object arg) {
        tv.setforeground(getSelectedItem());
        return true;
}
}

class PointChoice extends Choice {

TriangulationViewer tv;

PointChoice(TriangulationViewer tv) {
        super();
        this.tv = tv;
        addItem("3");
        addItem("4");
        addItem("5");
        addItem("10");
        addItem("20");
        addItem("50");
        addItem("100");
        addItem("500");
        addItem("1000");
        select("20");
        tv.numpoints = 20;
}

public boolean action (Event evt, Object arg) {
        tv.numpoints = Integer.valueOf((String)arg).intValue();
        tv.newrandom();
        return true;
}
}

class TriangulationCanvas extends Canvas {
        static final int POINTS = 1;
        static final int TRIANGLES = 2;
        static final int CIRCLES = 3;
        static final int CENTERS = 4;
        static final int BORDER = 5;
        static final int TESSELLATION = 6;

        static final int NODRAW = 0;
        static final int DRAWOVAL = 1;
        static final int FILLOVAL = 2;
        static final int DRAWTRIANGLE = 3;
        static final int FILLTRIANGLE = 4;

        Applet applet;
        Delaunay D;
        int     points;
        int     triangles;
        int     display;
        double  xmin, xmax;
        double  ymin, ymax;
        double  range, xoff, yoff;
        boolean paintall = true;
        int     highlight = NODRAW;

        public void setTriangulation(Delaunay d) {
                D = d;
                localupdate();
                paintall = true;
                repaint();
        }
        public void setDisplay(int i) {
                display = i;
                paintall = true;
                repaint();
        }
        public void setforeground(int i) {
                highlight = i;
                paintall = true;
                repaint();
        }
        public TriangulationCanvas (Applet a, Delaunay d) {
                this.applet = a;
                setDisplay(POINTS);
                setTriangulation(d);
                paintall = true;
                localupdate();
        }
        private void localupdate () {
                points = D.getNumberPoints();
                triangles = D.getNumberTriangles();
                xmin = D.getXMin();
                xmax = D.getXMax();
                ymin = D.getYMin();
                ymax = D.getYMax();
                double xdel = xmax - xmin;
                double ydel = ymax - ymin;
                range = Math.max(xdel,ydel) * 1.2;
                xoff = (xmin + xmax - range) / 2;
                yoff = (ymin + ymax - range) / 2;
        }

        private double[] center(double x1, double y1,
                                double x2, double y2,
                                double x3, double y3) {
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

        private int ccw (double p0[], double p1[], double p2[]) {
                double dx1, dx2, dy1, dy2, d;
                dx1 = p1[0] - p0[0];
                dy1 = p1[1] - p1[1];
                dx2 = p2[0] - p0[0];
                dy2 = p2[1] - p0[1];
                d = dx1*dy2 - dy1*dx2;
                if (d > 0) return 1;
                if (d < 0) return -1;
                if (dx1*dx2<0 || dy1*dy2 < 0) return -1;
                if (dx1*dx1+dy1*dy1 < dx2*dx2+dy2*dy2) return 1;
                return 0;
        }

        private boolean onright(double p[], double p1[], double p2[]) {
                return (p1[1] - p[1]) * (p2[0] - p[0]) >
                       (p1[0] - p[0]) * (p2[1] - p[1]);
        }

        private double[] center(int v[]) {
                double p1[] = D.getPoint(v[0]);
                double p2[] = D.getPoint(v[1]);
                double p3[] = D.getPoint(v[2]);
                double c[] = center(p1[0],p1[1],p2[0],p2[1],
                                    p3[0],p3[1]);
                return c;
        }

        private double[] center (int i) {
                return center(D.getTriangle(i));
        }

        private double[] midpoint (double p1[], double p2[]) {
                double m[] = new double[2];
                m[0] = (p1[0] + p2[0]) / 2;
                m[1] = (p1[1] + p2[1]) / 2;
                return m;
        }

        private void longLine(Graphics g, double c[], double m[], boolean dir) {
                int nx1, ny1, nx2, ny2;
                double x1 = wo + (c[0]-xoff)*r;
                double y1 = ho - (c[1]-yoff)*r;
                double x2 = wo + (m[0]-xoff)*r;
                double y2 = ho - (m[1]-yoff)*r;

                if (x1 <= 0 || x1 >= width ||
                    y1 <= 0 || y1 >= height)
                        return;
                dir = !dir;
                if (x1 == x2) {
                        if (y1 < y2 == dir) {
                                nx1 = (int)x1;
                                ny1 = (int)y1;
                                nx2 = (int)x2;
                                ny2 = height;
                        } else {
                                nx1 = (int)x1;
                                ny1 = (int)y1;
                                nx2 = (int)x2;
                                ny2 = 0;
                        }
                } else if (x1 < x2 == dir) {
                                nx1 = (int)x1;
                                ny1 = (int)y1;
                                nx2 = (int)width;
                                ny2 = (int)((y2-y1)*(width-x1)/(x2-x1)+y1);
                } else {
                                nx1 = (int)x1;
                                ny1 = (int)y1;
                                nx2 = (int)0;
                                ny2 = (int)((y2-y1)*(0-x1)/(x2-x1)+y1);
                }
                g.drawLine(nx1,ny1,nx2,ny2);
        }

        int     lastt = -1, newt = -1;
        public void update (Graphics g) {
                if (paintall == true) {
                        paint(g);
                        newt = -1;
                } else if (lastt != newt) {
                        foreground(g,lastt);
                        foreground(g,newt);
                }
                lastt = newt;
        }

        private void foreground (Graphics g, int tri) {
                        if (highlight == NODRAW) return;
                        if (tri < 0) return;
                        int v[] = D.getTriangle(tri);
                        double p2[];
                        g.setXORMode(Color.gray);
                        switch (highlight) {
                        case DRAWOVAL:
                        case FILLOVAL:
                                p2 = D.getPoint(v[1]);
                                double c[] = center(v);
                                int cx = wo + (int)((c[0]-xoff)*r);
                                int cy = ho - (int)((c[1]-yoff)*r);
                                double dx = c[0] - p2[0];
                                double dy = c[1] - p2[1];
                                double rad = Math.sqrt(dx*dx+dy*dy);
                                int ir = (int)(rad*r);
                                if (highlight == DRAWOVAL) {
                                        g.fillOval(cx-2,cy-2,5,5);
                                        g.drawOval(cx-ir,cy-ir,2*ir,2*ir);
                                } else {
                                        g.fillOval(cx-ir,cy-ir,2*ir,2*ir);
                                }
                                break;
                        case DRAWTRIANGLE:
                        case FILLTRIANGLE:
                                double p1[] = D.getPoint(v[0]);
                                p2 = D.getPoint(v[1]);
                                double p3[] = D.getPoint(v[2]);
                                int x[] = new int[4];
                                int y[] = new int[4];
                                x[0] = wo + (int)((p1[0]-xoff)*r);
                                y[0] = ho - (int)((p1[1]-yoff)*r);
                                x[1] = wo + (int)((p2[0]-xoff)*r);
                                y[1] = ho - (int)((p2[1]-yoff)*r);
                                x[2] = wo + (int)((p3[0]-xoff)*r);
                                y[2] = ho - (int)((p3[1]-yoff)*r);
                                x[3] = x[0];
                                y[3] = y[0];
                                if (highlight == DRAWTRIANGLE) {
                                        g.drawPolygon(x,y,4);
                                } else {
                                        g.fillPolygon(x,y,3);
                                }
                                break;
                        }
                        g.setPaintMode();
        }

        int wo, ho;
        double r;
        int     width, height;

        public void paint(Graphics g) {
                Dimension s = size();
                width = s.width;
                height = s.height;
                int     r1 = Math.min(width, height);
                wo = (width - r1) / 2;
                ho = height - (height - r1) / 2;
                r = r1 / range;

                g.setColor(getBackground());
//                g.fill3DRect(0, 0, width, height, true);
                g.fillRect(0, 0, width, height);
                switch (display) {
                case POINTS:
                        g.setColor(Color.blue);
                        for (int i = 0; i<points;i++) {
                                double p[] = D.getPoint(i);
                                int x = wo + (int)((p[0]-xoff)*r);
                                int y = ho - (int)((p[1]-yoff)*r);
                                g.drawString(""+i,x+1,y);
                        }
                        break;
                case TRIANGLES:
                        g.setColor(Color.black);
                        for (int i = 0; i<triangles;i++) {
                                int v[] = D.getTriangle(i);
                                double c[] = center(v);
                                double p1[] = D.getPoint(v[0]);
                                double p2[] = D.getPoint(v[1]);
                                double p3[] = D.getPoint(v[2]);
                                int x1 = wo + (int)((p1[0]-xoff)*r);
                                int y1 = ho - (int)((p1[1]-yoff)*r);
                                int x2 = wo + (int)((p2[0]-xoff)*r);
                                int y2 = ho - (int)((p2[1]-yoff)*r);
                                int x3 = wo + (int)((p3[0]-xoff)*r);
                                int y3 = ho - (int)((p3[1]-yoff)*r);
                                g.drawLine(x1,y1,x2,y2);
                                g.drawLine(x2,y2,x3,y3);
                                g.drawLine(x3,y3,x1,y1);
                        }
                        break;
                case CIRCLES:
                        g.setColor(Color.yellow);
                        for (int i = 0; i<triangles;i++) {
                                int v[] = D.getTriangle(i);
                                double p2[] = D.getPoint(v[1]);
                                double c[] = center(v);
                                int cx = wo + (int)((c[0]-xoff)*r);
                                int cy = ho - (int)((c[1]-yoff)*r);
                                g.drawLine(cx,cy,cx,cy);
                                double dx = c[0] - p2[0];
                                double dy = c[1] - p2[1];
                                double rad = Math.sqrt(dx*dx+dy*dy);
                                int ir = (int)(rad*r);
                                g.drawOval(cx-ir,cy-ir,2*ir,2*ir);
                        }
                        break;
                case CENTERS:
                        g.setColor(Color.black);
                        for (int i = 0; i<triangles;i++) {
                                double c[] = center(i);
                                int cx = wo + (int)((c[0]-xoff)*r);
                                int cy = ho - (int)((c[1]-yoff)*r);
                                g.drawLine(cx,cy,cx,cy);
                        }
                        break;
                case BORDER:
                        g.setColor(Color.black);
                        for (int i = 0; i<triangles;i++) {
                                int v[] = D.getTriangle(i);
                                double p1[] = D.getPoint(v[0]);
                                double p2[] = D.getPoint(v[1]);
                                double p3[] = D.getPoint(v[2]);
                                int x1 = wo + (int)((p1[0]-xoff)*r);
                                int y1 = ho - (int)((p1[1]-yoff)*r);
                                int x2 = wo + (int)((p2[0]-xoff)*r);
                                int y2 = ho - (int)((p2[1]-yoff)*r);
                                int x3 = wo + (int)((p3[0]-xoff)*r);
                                int y3 = ho - (int)((p3[1]-yoff)*r);

                                int e[] = D.getEdge(i);
                                if (e[0] == -1)
                                        g.drawLine(x1,y1,x2,y2);
                                if (e[1] == -1)
                                        g.drawLine(x2,y2,x3,y3);
                                if (e[2] == -1)
                                        g.drawLine(x3,y3,x1,y1);
                        }
                        break;
                case TESSELLATION:
                        g.setColor(Color.black);
                        for (int i = 0; i<triangles;i++) {
                                int v[] = D.getTriangle(i);
                                double p1[] = D.getPoint(v[0]);
                                double p2[] = D.getPoint(v[1]);
                                double p3[] = D.getPoint(v[2]);
                                int e[] = D.getEdge(i);
                                double c1[] = center(v);
                                boolean cr;
                                int c1x = wo + (int)((c1[0]-xoff)*r);
                                int c1y = ho - (int)((c1[1]-yoff)*r);

                                if (e[0] == -1) {
                                        double m[] = midpoint(p1,p2);
                                        cr = onright(c1, p1, p2);
                                        longLine(g,c1,m,cr);
                                }
                                if (e[1] == -1) {
                                        double m[] = midpoint(p2,p3);
                                        cr = onright(c1, p2, p3);
                                        longLine(g,c1,m,cr);
                                }
                                if (e[2] == -1) {
                                        double m[] = midpoint(p3,p1);
                                        cr = onright(c1, p3, p1);
                                        longLine(g,c1,m,cr);
                                }
                                for (int j = 0; j < 3; j++) {
                                        if (e[j] == -1) continue;
                                        double c2[] = center(e[j]);
                                        int c2x = wo + (int)((c2[0]-xoff)*r);
                                        int c2y = ho - (int)((c2[1]-yoff)*r);
                                        g.drawLine(c1x,c1y,c2x,c2y);
                                }
                        }
                        break;
                default: // unknown
                        break;
                }

                // draw centers of circles
                g.setColor(Color.yellow);
                for (int i = 0; i<triangles;i++) {
                        int v[] = D.getTriangle(i);
                        double c1[] = center(v);
                        int c1x = wo + (int)((c1[0]-xoff)*r);
                        int c1y = ho - (int)((c1[1]-yoff)*r);
                        g.fillOval(c1x-2,c1y-2,5,5);
                }

                // finally, draw points in the data set
                g.setColor(Color.red);
                for (int i = 0; i<points;i++) {
                        double p[] = D.getPoint(i);
                        int x = wo + (int)((p[0]-xoff)*r);
                        int y = ho - (int)((p[1]-yoff)*r);
                        g.fillOval(x-2,y-2,5,5);
                }
                paintall = false;
        }

        private void handlepoint(int x, int y) {
                Dimension s = size();
                int     width = s.width, height = s.height;
                int     r1 = Math.min(width, height);
                int     wo = (width - r1) / 2;
                int     ho = height - (height - r1) / 2;
                double  r = r1 / range;
                double  px, py;

                px = (x - wo)  / r + xoff;
                py = (-y + ho) / r + yoff;
                try {
                        newt = D.locate(px, py);
                } catch (Exception e) { newt = -1; }
                if (newt >= 0) {
                        applet.showStatus("("+px+","+py+") in triangle "+newt);
                } else {
                        applet.showStatus("("+px+","+py+")");
                }
                repaint();
        }

        public boolean handleEvent(Event evt) {
                if (evt.id == Event.MOUSE_MOVE) {
                        handlepoint(evt.x, evt.y);
                        return true;
                }
                return false;
        }
}
