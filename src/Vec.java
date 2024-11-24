import java.text.DecimalFormat;

public class Vec {
    double x = 0, y = 0, z = 0;

    public Vec(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public Vec() {
    }

    public Vec addVec(Vec vec) {
        return new Vec(x + vec.x, y + vec.y, z + vec.z);
    }

    public Vec subVec(Vec vec) {
        return new Vec(x - vec.x, y - vec.y, z - vec.z);
    }

    public Vec mul(double m) {
        return new Vec(x * m, y * m, z * m);
    }

    public Vec mul(Matrix mtr){return new Vec();}

    public Vec div(double d) {
        if (d == 0){
            throw new IllegalArgumentException("division by zero");
        }
        return new Vec(x / d, y / d, z / d);
    }

    double distance() {
        return Math.sqrt(x * x + y * y + z * z);
    }

    public Vec(Vec vec) {
        this.x = vec.x;
        this.y = vec.y;
        this.z = vec.z;
    }

    @Override
    public String toString() {
        return "{" +
                x +
                ", " + y +
                ", " + z +
                '}';
    }

    public String toShortString(){
        DecimalFormat df = new DecimalFormat("#.##");
        return "{" +
                df.format(x) +
                ", " + df.format(y) +
                ", " + df.format(z) +
                '}';
    }
}
