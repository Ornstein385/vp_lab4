import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class Main {

    static double a = 0.2, b = 0.2, c = 5.7, minValueHor = 0.0, maxValueHor = 0.5, minValueVer = Double.MAX_VALUE,
            maxValueVer = Double.MIN_VALUE, h = 0.1, q = 0.5, yAxisBegin = 0, yAxisInterval = 2.5;

    static char selectedVariable = 'x', selectedParameter = 'a';

    static int resolution = 10000, stepCount = 10000;

    static double[][] diagramValues;

    static JFrame frame = new JFrame("КТПрВП лабораторная");

    static JComboBox<String> methodSelection = new JComboBox<>(new String[]{"Полуявный метод Эйлера-Кромера", "Явный метод Эйлера",
            "Полунеявный метод CD", "Неявный метод Эйлера", "Явный метод Рунге-Кутты RK4", "Явный метод Рунге-Кутты RK3/8",
            "Композиционная схема s3ord4", "Композиционная схема s7ord6", "Композиционная схема s9ord8"
    }),
            problemSelection = new JComboBox<>(new String[]{"Аттрактор Рёсслера", "Аттрактор Нозе-Гувера"}),
            parameterSelection = new JComboBox<>(new String[]{"a", "b", "c"}),
            variableSelection = new JComboBox<>(new String[]{"x", "y", "z"});

    static JLabel parametersLabel = new JLabel("параметры системы (a, b, c):"),
            intervalLabel = new JLabel("строим график с параметром                    в интервале от                        до                  по переменной"),
            label1 = new JLabel("разрешение"),
            label2 = new JLabel("количество шагов интегрирования"),
            label3 = new JLabel("шаг интегрирования"),
            label4 = new JLabel("считаем режим установившемся после __% шагов"),
            label5 = new JLabel("начальная координата по вертикали"),
            label6 = new JLabel("длина единицы шкалы по вертикали"),
            applied = new JLabel("");

    static JTextField fromValueField = new JTextField("0"), toValueField = new JTextField("0.5"),
            aField = new JTextField("0.2"), bField = new JTextField("0.2"), cField = new JTextField("5.7"),
            field1 = new JTextField("10000"),
            field2 = new JTextField("10000"),
            field3 = new JTextField("0.1"),
            field4 = new JTextField("50"),
            field5 = new JTextField("0"),
            field6 = new JTextField("2.5");

    static JButton apply = new JButton("применить");

    static ArrayList<Double> getPeaks(ArrayList<Vec> values) {
        ArrayList<Double> peaks = new ArrayList<>();
        for (int i = 1; i < values.size() - 1; i++) {
            switch (selectedVariable) {
                case 'x':
                    if (values.get(i - 1).x < values.get(i).x && values.get(i + 1).x < values.get(i).x) {
                        peaks.add(values.get(i).x);
                    }
                    break;
                case 'y':
                    if (values.get(i - 1).y < values.get(i).y && values.get(i + 1).y < values.get(i).y) {
                        peaks.add(values.get(i).y);
                    }
                    break;
                case 'z':
                    if (values.get(i - 1).z < values.get(i).z && values.get(i + 1).z < values.get(i).z) {
                        peaks.add(values.get(i).z);
                    }
                    break;
            }

        }
        return peaks;
    }

    static ArrayList<Vec> implicitEulerMethod_Rossler() {
        ArrayList<Vec> values = new ArrayList<>();
        Vec vars = new Vec();
        for (int i = 0; i < stepCount; i++) {
            Vec m = new Vec(vars);
            for (int j = 0; j < 10 || m.distance() > 1e-14; j++) {
                Matrix e = new Matrix(new double[][]{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.1}});
                Matrix jackobian = new Matrix(new double[][]{{0, -1, -1}, {1, a, 0}, {vars.z, 0, vars.x - c}});
                Matrix inv = e.sub(jackobian).inversion();
                m = m.subVec(vars).subVec(new Vec(-vars.y - vars.z, vars.x + a * vars.y, b + vars.z * (vars.x - c))).mul(inv);
            }
            vars.addVec(m);
            if (!Double.isNaN(vars.x) && !Double.isNaN(vars.y) && !Double.isNaN(vars.z)) {
                values.add(vars);
            } else {
                break;
            }
        }
        return values;
    }

    static ArrayList<Vec> implicitEulerMethod_NoseHoover() {
        ArrayList<Vec> values = new ArrayList<>();
        Vec vars = new Vec();
        for (int i = 0; i < stepCount; i++) {
            Vec m = new Vec(vars);
            for (int j = 0; j < 10 || m.distance() > 1e-14; j++) {
                Matrix e = new Matrix(new double[][]{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.1}});
                Matrix jackobian = new Matrix(new double[][]{{0, a, 0}, {-1, vars.z, vars.y}, {0, -2 * vars.y, 0}});
                Matrix inv = e.sub(jackobian).inversion();
                m = m.subVec(vars).subVec(new Vec(a * vars.y, -vars.x + vars.y * vars.z, 1 - vars.y * vars.y)).mul(inv);
            }
            vars.addVec(m);
            if (!Double.isNaN(vars.x) && !Double.isNaN(vars.y) && !Double.isNaN(vars.z)) {
                values.add(vars);
            } else {
                break;
            }
        }
        return values;
    }

    static Vec nextStepRossler(Vec vec) {
        double x = vec.x, y = vec.y, z = vec.z;
        x += h / 2 * (-y - z);
        y += h / 2 * (x + a * y);
        z += h / 2 * (b + z * (x - c));
        z += h / 2 * (b + z * (x - c));
        y += h / 2 * (x + a * y);
        x += h / 2 * (-y - z);
        return new Vec(x - vec.x, y - vec.y, z - vec.z);
    }

    static Vec nextStepNoseHoover(Vec vec) {
        double x = vec.x, y = vec.y, z = vec.z;
        x += h / 2 * (a * y);
        y += h / 2 * (-x + y * z);
        z += h / 2 * (1 - y * y);
        z += h / 2 * (1 - y * y);
        y += h / 2 * (-x + y * z);
        x += h / 2 * (a * y);
        return new Vec(x - vec.x, y - vec.y, z - vec.z);
    }

    static ArrayList<Vec> rk4_Rossler() {
        ArrayList<Vec> values = new ArrayList<>();
        Vec vars = new Vec();
        for (int i = 0; i < stepCount; i++) {
            Vec[] k = new Vec[4];
            k[0] = nextStepRossler(vars);
            k[1] = nextStepRossler(k[0].mul(0.5).addVec(vars));
            k[2] = nextStepRossler(k[1].mul(0.5).addVec(vars));
            k[3] = nextStepRossler(k[2].addVec(vars));
            vars = vars.addVec(k[0].mul(1.0/6)).addVec(k[1].mul(1.0/3)).addVec(k[2].mul(1.0/3)).addVec(k[3].mul(1.0/6));
            if (i > stepCount * q) {
                values.add(vars);
            }
        }
        return values;
    }

    static ArrayList<Vec> rk4_NoseHoover() {
        ArrayList<Vec> values = new ArrayList<>();
        Vec vars = new Vec(1e-6, 1e-6, 1e-6);
        for (int i = 0; i < stepCount; i++) {
            Vec[] k = new Vec[4];
            k[0] = nextStepNoseHoover(vars);
            k[1] = nextStepNoseHoover(k[0].mul(0.5).addVec(vars));
            k[2] = nextStepNoseHoover(k[1].mul(0.5).addVec(vars));
            k[3] = nextStepNoseHoover(k[2].addVec(vars));
            vars = vars.addVec(k[0].mul(1.0/6)).addVec(k[1].mul(1.0/3)).addVec(k[2].mul(1.0/3)).addVec(k[3].mul(1.0/6));
            if (i > stepCount * q) {
                values.add(vars);
            }
        }
        return values;
    }

    static ArrayList<Vec> rk38_Rossler() {
        ArrayList<Vec> values = new ArrayList<>();
        Vec vars = new Vec();
        for (int i = 0; i < stepCount; i++) {
            Vec[] k = new Vec[4];
            k[0] = nextStepRossler(vars);
            k[1] = nextStepRossler((k[0].mul(1.0/3)).addVec(vars));
            k[2] = nextStepRossler((k[0].mul(-1.0/3)).addVec(k[1]).addVec(vars));
            k[3] = nextStepRossler(k[0].subVec(k[1]).addVec(k[2]).addVec(vars));
            vars = vars.addVec(k[0].mul(1.0/8)).addVec(k[1].mul(3.0/8)).addVec(k[2].mul(3.0/8)).addVec(k[3].mul(1.0/8));
            if (i > stepCount * q) {
                values.add(vars);
            }
        }
        return values;
    }

    static ArrayList<Vec> rk38_NoseHoover() {
        ArrayList<Vec> values = new ArrayList<>();
        Vec vars = new Vec(1e-6, 1e-6, 1e-6);
        for (int i = 0; i < stepCount; i++) {
            Vec[] k = new Vec[4];
            k[0] = nextStepNoseHoover(vars);
            k[1] = nextStepNoseHoover((k[0].mul(1.0/3)).addVec(vars));
            k[2] = nextStepNoseHoover((k[0].mul(-1.0/3)).addVec(k[1]).addVec(vars));
            k[3] = nextStepNoseHoover(k[0].subVec(k[1]).addVec(k[2]).addVec(vars));
            vars = vars.addVec(k[0].mul(1.0/8)).addVec(k[1].mul(3.0/8)).addVec(k[2].mul(3.0/8)).addVec(k[3].mul(1.0/8));
            if (i > stepCount * q) {
                values.add(vars);
            }
        }
        return values;
    }




    static ArrayList<Vec> comp_s3ord4_Rossler() {
        ArrayList<Vec> values = new ArrayList<>();
        double[] k = new double[]{1.3512071919596578, -1.7024143839193155, 1.3512071919596578};
        double x = 0, y = 0, z = 0;
        for (int i = 0; i < stepCount; i++) {
            for (int j = 0; j < k.length; j++){
                x += k[j] * h / 2 * (-y - z);
                y += k[j] * h / 2 * (x + a * y);
                z += k[j] * h / 2 * (b + z * (x - c));
                z += k[j] * h / 2 * (b + z * (x - c));
                y += k[j] * h / 2 * (x + a * y);
                x += k[j] * h / 2 * (-y - z);
            }
            if (i > stepCount * q && !Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> comp_s7ord6_Rossler() {
        ArrayList<Vec> values = new ArrayList<>();
        double[] k = new double[]{0.78451361047755726382,0.23557321335935813368,-1.1776799841788710069,
                1.3151863206839112189, -1.1776799841788710069,0.23557321335935813368,0.78451361047755726382};
        double x = 0, y = 0, z = 0;
        for (int i = 0; i < stepCount; i++) {
            for (int j = 0; j < k.length; j++){
                x += k[j] * h / 2 * (-y - z);
                y += k[j] * h / 2 * (x + a * y);
                z += k[j] * h / 2 * (b + z * (x - c));
                z += k[j] * h / 2 * (b + z * (x - c));
                y += k[j] * h / 2 * (x + a * y);
                x += k[j] * h / 2 * (-y - z);
            }
            if (i > stepCount * q && !Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> comp_s9ord8_Rossler() {
        ArrayList<Vec> values = new ArrayList<>();
        double[] k = new double[]{0.39103020330868478817,0.33403728961113601749,-0.70622728118756134346,0.08187754964805944576890,
                0.79856447723936218406,0.08187754964805944576890,-0.70622728118756134346,0.33403728961113601749,0.39103020330868478817};
        double x = 0, y = 0, z = 0;
        for (int i = 0; i < stepCount; i++) {
            for (int j = 0; j < k.length; j++){
                x += k[j] * h / 2 * (-y - z);
                y += k[j] * h / 2 * (x + a * y);
                z += k[j] * h / 2 * (b + z * (x - c));
                z += k[j] * h / 2 * (b + z * (x - c));
                y += k[j] * h / 2 * (x + a * y);
                x += k[j] * h / 2 * (-y - z);
            }
            if (i > stepCount * q && !Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> comp_s3ord4_NoseHoover() {
        ArrayList<Vec> values = new ArrayList<>();
        double[] k = new double[]{1.3512071919596578, -1.7024143839193155, 1.3512071919596578};
        double x = 1e-6, y = 1e-6, z = 1e-6;
        for (int i = 0; i < stepCount; i++) {
            for (int j = 0; j < k.length; j++){
                x += k[j] * h / 2 * (a * y);
                y += k[j] * h / 2 * (-x + y * z);
                z += k[j] * h / 2 * (1 - y * y);
                z += k[j] * h / 2 * (1 - y * y);
                y += k[j] * h / 2 * (-x + y * z);
                x += k[j] * h / 2 * (a * y);
            }
            if (i > stepCount * q && !Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> comp_s7ord6_NoseHoover() {
        ArrayList<Vec> values = new ArrayList<>();
        double[] k = new double[]{0.78451361047755726382,0.23557321335935813368,-1.1776799841788710069,
                1.3151863206839112189, -1.1776799841788710069,0.23557321335935813368,0.78451361047755726382};
        double x = 1e-6, y = 1e-6, z = 1e-6;
        for (int i = 0; i < stepCount; i++) {
            for (int j = 0; j < k.length; j++){
                x += k[j] * h / 2 * (a * y);
                y += k[j] * h / 2 * (-x + y * z);
                z += k[j] * h / 2 * (1 - y * y);
                z += k[j] * h / 2 * (1 - y * y);
                y += k[j] * h / 2 * (-x + y * z);
                x += k[j] * h / 2 * (a * y);
            }
            if (i > stepCount * q && !Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> comp_s9ord8_NoseHoover() {
        ArrayList<Vec> values = new ArrayList<>();
        double[] k = new double[]{0.39103020330868478817,0.33403728961113601749,-0.70622728118756134346,0.08187754964805944576890,
                0.79856447723936218406,0.08187754964805944576890,-0.70622728118756134346,0.33403728961113601749,0.39103020330868478817};
        double x = 1e-6, y = 1e-6, z = 1e-6;
        for (int i = 0; i < stepCount; i++) {
            for (int j = 0; j < k.length; j++){
                x += k[j] * h / 2 * (a * y);
                y += k[j] * h / 2 * (-x + y * z);
                z += k[j] * h / 2 * (1 - y * y);
                z += k[j] * h / 2 * (1 - y * y);
                y += k[j] * h / 2 * (-x + y * z);
                x += k[j] * h / 2 * (a * y);
            }
            if (i > stepCount * q && !Double.isNaN(x) && !Double.isNaN(y) && !Double.isNaN(z)) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }





    static ArrayList<Vec> explicitEulerMethod_Rossler() {
        ArrayList<Vec> values = new ArrayList<>();
        double x = 0, y = 0, z = 0;
        for (int i = 0; i < stepCount; i++) {
            double dx = h * (-y - z);
            double dy = h * (x + a * y);
            double dz = h * (b + z * (x - c));
            x += dx;
            y += dy;
            z += dz;
            if (i > stepCount * q) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> semiExplicitEulerMethod_Rossler() {
        ArrayList<Vec> values = new ArrayList<>();
        double x = 0, y = 0, z = 0;
        for (int i = 0; i < stepCount; i++) {
            x += h * (-y - z);
            y += h * (x + a * y);
            z += h * (b + z * (x - c));
            if (i > stepCount * q) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> semiImplicitCDMethod_Rossler() {
        ArrayList<Vec> values = new ArrayList<>();
        double x = 0, y = 0, z = 0;
        for (int i = 0; i < stepCount; i++) {
            x += h / 2 * (-y - z);
            y += h / 2 * (x + a * y);
            z += h / 2 * (b + z * (x - c));
            z += h / 2 * (b + z * (x - c));
            y += h / 2 * (x + a * y);
            x += h / 2 * (-y - z);
            if (i > stepCount * q) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> explicitEulerMethod_NoseHoover() {
        ArrayList<Vec> values = new ArrayList<>();
        double x = 1e-6, y = 1e-6, z = 1e-6;
        for (int i = 0; i < stepCount; i++) {
            double dx = h * (a * y);
            double dy = h * (-x + y * z);
            double dz = h * (1 - y * y);
            x += dx;
            y += dy;
            z += dz;
            if (i > stepCount * q) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> semiExplicitEulerMethod_NoseHoover() {
        ArrayList<Vec> values = new ArrayList<>();
        double x = 1e-6, y = 1e-6, z = 1e-6;
        for (int i = 0; i < stepCount; i++) {
            x += h * (a * y);
            y += h * (-x + y * z);
            z += h * (1 - y * y);
            if (i > stepCount * q) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> semiImplicitCDMethod_NoseHoover() {
        ArrayList<Vec> values = new ArrayList<>();
        double x = 1e-6, y = 1e-6, z = 1e-6;
        for (int i = 0; i < stepCount; i++) {
            x += h / 2 * (a * y);
            y += h / 2 * (-x + y * z);
            z += h / 2 * (1 - y * y);
            z += h / 2 * (1 - y * y);
            y += h / 2 * (-x + y * z);
            x += h / 2 * (a * y);
            if (i > stepCount * q) {
                values.add(new Vec(x, y, z));
            }
        }
        return values;
    }

    static ArrayList<Vec> implicitEulerMethod_Rosler() {
        return comp_s9ord8_Rossler();
    }

    static ArrayList<Vec> implicitEulerMethod_NoseHover() {
        return comp_s9ord8_NoseHoover();
    }

    public static void main(String[] args) {
        frame.setBounds(0, 0, 1920, 1080);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLayout(new BorderLayout(1, 1));
        frame.setVisible(true);

        apply.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    double updateA = Double.parseDouble(aField.getText());
                    double updateB = Double.parseDouble(bField.getText());
                    double updateC = Double.parseDouble(cField.getText());
                    double updateLeft = Double.parseDouble(fromValueField.getText());
                    double updateRight = Double.parseDouble(toValueField.getText());
                    int updateRes = Integer.parseInt(field1.getText());
                    int updateSteps = Integer.parseInt(field2.getText());
                    double updateH = Double.parseDouble(field3.getText());
                    double updateQ = Double.parseDouble(field4.getText()) / 100;
                    double updateYBegin = Double.parseDouble(field5.getText());
                    double updateInterval = Double.parseDouble(field6.getText());
                    if (updateLeft > updateRight || updateRes < 0 || updateSteps < 0 || updateH < 0 || updateQ < 0 || updateQ > 1 || updateInterval <= 0) {
                        throw new Exception();
                    }
                    a = updateA;
                    b = updateB;
                    c = updateC;
                    minValueHor = updateLeft;
                    maxValueHor = updateRight;
                    resolution = updateRes;
                    stepCount = updateSteps;
                    h = updateH;
                    q = updateQ;
                    yAxisBegin = updateYBegin;
                    yAxisInterval = updateInterval;
                    switch (parameterSelection.getSelectedIndex()) {
                        case 0 -> selectedParameter = 'a';
                        case 1 -> selectedParameter = 'b';
                        case 2 -> selectedParameter = 'c';
                    }
                    switch (variableSelection.getSelectedIndex()) {
                        case 0 -> selectedVariable = 'x';
                        case 1 -> selectedVariable = 'y';
                        case 2 -> selectedVariable = 'z';
                    }
                    applied.setText("параметры применены");
                    drawDiagram();
                } catch (Exception ex) {
                    applied.setText("параметры некорректны");
                    frame.repaint();
                }
            }
        });

        methodSelection.setBounds(0, 0, 230, 20);
        problemSelection.setBounds(0, 20, 230, 20);
        parametersLabel.setBounds(0, 40, 200, 20);
        aField.setBounds(0, 60, 50, 20);
        bField.setBounds(50, 60, 50, 20);
        cField.setBounds(100, 60, 50, 20);
        intervalLabel.setBounds(0, 80, 700, 20);
        parameterSelection.setBounds(180, 80, 40, 20);
        variableSelection.setBounds(550, 80, 40, 20);
        fromValueField.setBounds(340, 80, 40, 20);
        toValueField.setBounds(410, 80, 40, 20);
        label1.setBounds(0, 100, 500, 20);
        label2.setBounds(0, 120, 500, 20);
        label3.setBounds(0, 140, 500, 20);
        label4.setBounds(0, 160, 500, 20);
        field1.setBounds(300, 100, 70, 20);
        field2.setBounds(300, 120, 70, 20);
        field3.setBounds(300, 140, 70, 20);
        field4.setBounds(300, 160, 70, 20);
        apply.setBounds(0, 180, 100, 20);
        applied.setBounds(0, 200, 200, 20);
        label5.setBounds(250, 0, 300, 20);
        label6.setBounds(250, 20, 300, 20);
        field5.setBounds(500, 0, 50, 20);
        field6.setBounds(500, 20, 50, 20);

        frame.add(methodSelection);
        frame.add(problemSelection);
        frame.add(parametersLabel);
        frame.add(aField);
        frame.add(bField);
        frame.add(cField);
        frame.add(intervalLabel);
        frame.add(parameterSelection);
        frame.add(variableSelection);
        frame.add(fromValueField);
        frame.add(toValueField);
        frame.add(label1);
        frame.add(label2);
        frame.add(label3);
        frame.add(label4);
        frame.add(field1);
        frame.add(field2);
        frame.add(field3);
        frame.add(field4);
        frame.add(apply);
        frame.add(applied);
        frame.add(label5);
        frame.add(label6);
        frame.add(field5);
        frame.add(field6);

        frame.add(new Plane());
        drawDiagram();
    }

    static void drawDiagram() {
        diagramValues = new double[resolution][];
        switch (selectedParameter) {
            case 'a' -> a = minValueHor;
            case 'b' -> b = minValueHor;
            case 'c' -> c = minValueHor;
        }
        minValueVer = Double.MAX_VALUE;
        maxValueVer = Double.MIN_VALUE;
        for (int i = 0; i < resolution; i++) {
            ArrayList<Double> values = new ArrayList<>();
            switch (problemSelection.getSelectedIndex()) {
                case 0:
                    switch (methodSelection.getSelectedIndex()) {
                        case 0:
                            values = getPeaks(semiExplicitEulerMethod_Rossler());
                            break;
                        case 1:
                            values = getPeaks(explicitEulerMethod_Rossler());
                            break;
                        case 2:
                            values = getPeaks(semiImplicitCDMethod_Rossler());
                            break;
                        case 3:
                            values = getPeaks(implicitEulerMethod_Rosler());
                            break;
                        case 4:
                            values = getPeaks(rk4_Rossler());
                            break;
                        case 5:
                            values = getPeaks(rk38_Rossler());
                            break;
                        case 6:
                            values = getPeaks(comp_s3ord4_Rossler());
                            break;
                        case 7:
                            values = getPeaks(comp_s7ord6_Rossler());
                            break;
                        case 8:
                            values = getPeaks(comp_s9ord8_Rossler());
                            break;
                        default:
                            throw new RuntimeException();
                    }
                    break;
                case 1:
                    switch (methodSelection.getSelectedIndex()) {
                        case 0:
                            values = getPeaks(semiExplicitEulerMethod_NoseHoover());
                            break;
                        case 1:
                            values = getPeaks(explicitEulerMethod_NoseHoover());
                            break;
                        case 2:
                            values = getPeaks(semiImplicitCDMethod_NoseHoover());
                            break;
                        case 3:
                            values = getPeaks(implicitEulerMethod_NoseHover());
                            break;
                        case 4:
                            values = getPeaks(rk4_NoseHoover());
                            break;
                        case 5:
                            values = getPeaks(rk38_NoseHoover());
                            break;
                        case 6:
                            values = getPeaks(comp_s3ord4_NoseHoover());
                            break;
                        case 7:
                            values = getPeaks(comp_s7ord6_NoseHoover());
                            break;
                        case 8:
                            values = getPeaks(comp_s9ord8_NoseHoover());
                            break;
                        default:
                            throw new RuntimeException();
                    }
                    break;
            }
            diagramValues[i] = new double[values.size()];
            for (int j = 0; j < values.size(); j++) {
                diagramValues[i][j] = values.get(j);
                minValueVer = Math.min(diagramValues[i][j], minValueVer);
                maxValueVer = Math.max(diagramValues[i][j], maxValueVer);
            }
            switch (selectedParameter) {
                case 'a' -> a += (maxValueHor - minValueHor) / resolution;
                case 'b' -> b += (maxValueHor - minValueHor) / resolution;
                case 'c' -> c += (maxValueHor - minValueHor) / resolution;
            }
        }
        frame.repaint();
    }
}

class Plane extends JComponent {
    public void paintComponent(Graphics g) {
        DecimalFormat df = new DecimalFormat("#.###");
        //Main.maxValueVer = 50;
        //Main.minValueVer = 0;
        g.setColor(Color.BLACK);
        g.drawLine(100, 250, 100, 900);
        g.drawLine(100, 900, 1750, 900);
        for (int i = 0; i < 7; i++) {
            g.drawLine(100, 900 - i * 100, 80, 900 - i * 100);
            //g.drawString(df.format(Main.minValueVer + i * (Main.maxValueVer - Main.minValueVer) / 8), 50, 895 - i * 100);
            //g.drawString(df.format(i * 2.5), 50, 895 - i * 100);
            g.drawString(df.format(Main.yAxisBegin + i * Main.yAxisInterval), 50, 895 - i * 100);
        }
        for (int i = 0; i < 17; i++) {
            g.drawLine(100 + i * 100, 900, 100 + i * 100, 920);
            g.drawString(df.format(Main.minValueHor + i * (Main.maxValueHor - Main.minValueHor) / 16), 80 + i * 100, 935);
        }
        g.setColor(Color.BLUE);
        for (int i = 0; i < Main.diagramValues.length; i++) {
            for (int j = 0; j < Main.diagramValues[i].length; j++) {
                g.drawLine(100 + i * 1600 / Main.resolution, (int) (900 - (Main.diagramValues[i][j] - Main.yAxisBegin) * (100 / Main.yAxisInterval)), // (100 / ((Main.maxValueVer - Main.minValueVer) / 8))
                        100 + i * 1600 / Main.resolution, (int) (900 - (Main.diagramValues[i][j] - Main.yAxisBegin) * (100 / Main.yAxisInterval)));
            }
        }
    }
}