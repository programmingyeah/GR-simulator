using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
//using NumSharp;
using System.Diagnostics;
using System.Linq;

public class Tensor
{
    public int[] Shape;
    public double[] Data;

    public Tensor(params int[] shape)
    {
        int IndexCount = 1;

        this.Shape = shape;
        foreach (int dim in shape)
        {
            IndexCount *= dim;
        }

        Data = new double[IndexCount];


    }

    public Tensor(double[] data, params int[] shape)
    {
        Shape = shape;
        Data = data;

    }


    public double this[params int[] indeces]
    {
        get
        {
            if (indeces.Length == Shape.Length)
            {
                double returns;

                int funnyNumber = 0;
                int count = -1;
                foreach (int index in indeces)
                {

                    int funnyNumberTheSequel = 1;
                    for (int i = 0; i <= count; i++)
                    {

                        if (count != -1)
                        {
                            funnyNumberTheSequel *= Shape[count] - 1;
                        }
                        else
                        {
                            funnyNumberTheSequel = 1;
                        }

                    }

                    funnyNumber += index * funnyNumberTheSequel;
                    count++;
                }

                returns = Data[funnyNumber];

                return returns;
            }
            else
            {
                return 0;
            }
        }
        set
        {
            int funnyNumber = 0;
            int count = -1;
            foreach (int index in indeces)
            {
                
                int funnyNumberTheSequel = 1;
                for (int i = 0; i <= count; i++)
                {
                    
                    if (count != -1)
                    {
                        funnyNumberTheSequel *= Shape[count] - 1;
                    }
                    else
                    {
                        funnyNumberTheSequel = 1;
                    }
                    
                }

                funnyNumber += index * funnyNumberTheSequel;
                count++;
            }

            Data[funnyNumber] = value;
        }
    }

    public static Tensor operator +(Tensor t1, Tensor t2)
    {
        if (t1.Shape.SequenceEqual(t2.Shape) == false)
        {
            UnityEngine.Debug.Log("penis POOPOO EATER THE SHAPES DONT ALIGNNNN NIGGGAAA");
            return null;
        }
        else
        {
            Tensor t3 = new Tensor(t1.Shape);

            for (int i = 0; i < t1.Data.Length; i++)
            {
                t3.Data[i] = t1.Data[i] + t2.Data[i];
            }

            return t3;
        }
    }

    public static Tensor operator -(Tensor t1, Tensor t2)
    {
        if (t1.Shape.SequenceEqual(t2.Shape) == false)
        {
            UnityEngine.Debug.Log("penis POOPOO EATER THE SHAPES DONT ALIGNNNN NIGGGAAA");
            return null;
        }
        else
        {
            Tensor t3 = new Tensor(t1.Shape);

            for (int i = 0; i < t1.Data.Length; i++)
            {
                t3.Data[i] = t1.Data[i] - t2.Data[i];
            }

            return t3;
        }
    }

    public static Tensor operator *(double a, Tensor t)
    {
        Tensor T = new Tensor(t.Shape);

        for (int i = 0; i < t.Data.Length; i++)
        {
            T.Data[i] = a*t.Data[i];
        }

        return T;
    }

    public static Tensor operator *(Tensor t, double a)
    {
        Tensor T = new Tensor(t.Shape);

        for (int i = 0; i < t.Data.Length; i++)
        {
            T.Data[i] = a * t.Data[i];
        }

        return T;
    }

    public static Tensor operator *(Tensor t1, Tensor t2)
    {
        if (t1.Shape != t2.Shape)
        {
            UnityEngine.Debug.Log("you MORON THE SHSPES DONT ALIGNNNN");
            return null;
        }
        else
        {
            Tensor t3 = new Tensor(t1.Shape);

            for (int i = 0; i < t1.Data.Length; i++)
            {
                t3.Data[i] = t1.Data[i] * t2.Data[i];
            }

            return t3;
        }
    }

    public static Tensor operator /(Tensor t, double a)
    {

        Tensor T = new Tensor(t.Shape);

        for (int i = 0; i < t.Data.Length; i++)
        {
            T.Data[i] = t.Data[i]/a;
        }

        return T;
    }

    public static Tensor operator -(Tensor t)
    {
        Tensor T = new Tensor(t.Shape);

        for (int i = 0; i < t.Data.Length; i++)
        {
            T.Data[i] = -t.Data[i];
        }

        return T;
    }

    public void display()
    {
        string values = string.Join(", ", Data);
        UnityEngine.Debug.Log(values);
    }
}

public class GeodesicComputer : MonoBehaviour
{

    [SerializeField]
    Transform blackHole;

    [SerializeField]
    Transform body;


    const double G = 6.67408e-11;
    const double c = 299792458;
    const double scale = 1.496E6; //1.496E11

    Tensor schwarzschildMetric(Tensor coords, double mass)
    {

        Tensor g = new Tensor(3,3);

        double rs = 2 * G * mass / (c*c);
        double p = 1 - rs/coords[1];

        g[0, 0] = -c * c * p;
        g[1, 1] = 1 / p;
        g[2, 2] = coords[1] * coords[1];

        return g;
    }

    Tensor inverseMetric(Tensor coords, double mass)
    {

        Tensor g = new Tensor(3,3);

        double rs = 2 * G * mass / (c * c);
        double p = 1 - rs / coords[1];

        g[0, 0] = 1/(-c * c * p);
        g[1, 1] = p;
        g[2, 2] = 1/(coords[1] * coords[1]);

        return g;
    }

    Tensor christoffelSymbols(Tensor g, Tensor coords, double mass)
    {
        Stopwatch stopwatch = new Stopwatch();
        Tensor Gamma = new Tensor(3,3,3);
        Tensor inv = inverseMetric(coords, mass);
        Tensor[] metricDerivatives = {
            new Tensor(3,3),
            new Tensor(3,3),
            new Tensor(3,3),
        };
        double h = scale/10000.0;

        for (int i = 0; i < 3; i++)
        {
            Tensor dCoords = new Tensor(6);
            dCoords[i] = h;
            metricDerivatives[i] = (schwarzschildMetric(coords + dCoords, mass) - g) / h;
        }

        stopwatch.Start();
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                        for (int a = 0; a < 3; a++)
                        {

                            //double gDj = metricDerivatives[j, a, k];
                            //double gDk = metricDerivatives[k, a, j];
                            //double gDa = metricDerivatives[a, j, k];

                            Gamma[i, j, k] += 0.5 * inv[i, a] * (metricDerivatives[k][a, j] + metricDerivatives[j][a, k] - metricDerivatives[a][j, k]);

                        }
                }
            }
        }

        stopwatch.Stop();
        //UnityEngine.Debug.Log("Connection computation time: " + stopwatch.ElapsedMilliseconds + " ms");

        return Gamma;
    }

    Tensor cartToPol(Tensor cart)
    {
        double r = Math.Sqrt(cart[1] * cart[1] + cart[2] * cart[2]);

        Tensor result = cart;
        result[1] = r;
        double theta;
        if (cart[1] != 0)
        { 
            theta = Math.Atan(cart[2] / cart[1]);

        }
        else
        {
            if (cart[2] > 0)
            {
                theta = Math.PI / 2;
            }
            else
            {
                theta = 3 * Math.PI / 2;
            }
        }

        result[2] = theta;
        return result;
    }

    Tensor polToCart(Tensor pol)
    {
        Tensor result = new Tensor(3);

        result[0] = pol[0];
        result[1] = pol[1] * Math.Cos(pol[2]);
        result[2] = pol[1] * Math.Sin(pol[2]);

        return result;
    }

    Tensor F(Tensor position, double mass)
    {
        Stopwatch stopwatch = new Stopwatch();
        //stopwatch.Start();
        Tensor g = schwarzschildMetric(position, mass);
        //g.display();
        //stopwatch.Stop();
        //stopwatch.Start();
        Tensor Gamma = christoffelSymbols(g, position, mass);
        //UnityEngine.Debug.Log(Gamma[0,0,1]);
        //UnityEngine.Debug.Log(Gamma[1, 0, 0]);
        //UnityEngine.Debug.Log(Gamma[1, 2, 2]);
        //stopwatch.Stop();
        //UnityEngine.Debug.Log("Connection computation time: " + stopwatch.ElapsedMilliseconds + " ms");

        Tensor result = new Tensor(6);
        Tensor v = new Tensor(3);
        v[0] = position[3];
        v[1] = position[4];
        v[2] = position[5];
        result[0] = v[0];
        result[1] = v[1];
        result[2] = v[2];

        stopwatch.Start();
        for (int a = 0; a < 3; a++)
        {

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    result[a+3] -= Gamma[a, i, j] * v[i] * v[j];
                }
            }
        }
        stopwatch.Stop();

        //UnityEngine.Debug.Log("Connection computation time: " + stopwatch.ElapsedMilliseconds + " ms");
        //Debug.Log(result);
        return result;
    }

    Tensor solveGeodesic(double h, Tensor coords, double mass)
    {
        Tensor returnValues = new Tensor(6);
        //returnValues = coords + h * F(coords, mass);
        //Debug.Log(returnValues);
        Tensor k1 = F(coords, mass);
        Tensor k2 = F(coords + k1 * h / 2, mass);
        Tensor k3 = F(coords + k2 * h / 2, mass);
        Tensor k4 = F(coords + k3 * h, mass);

        returnValues = coords + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;

        return returnValues;
    }

    Tensor X = new Tensor(6);

    void Start()
    {
        float schwarzschildRadius = (float)(2*G*M/(c*c*scale));
        blackHole.localScale = new Vector3(schwarzschildRadius, schwarzschildRadius, 1);

        X[0] = 0;
        X[1] = scale * body.localPosition.x;
        X[2] = scale * body.localPosition.y;
        X[3] = 1;
        X[4] = 0;
        X[5] = 100; //1.99236428E-7
        X = cartToPol(X);

        //double h = 10000;

        //X = solveGeodesic(h, X, M);
        //X.display();
        // X = solveGeodesic(h, X, M);
        //X = solveGeodesic(h, X, M);
        //X = solveGeodesic(h, X, M);
        //X = solveGeodesic(h, X, M);
        //X = solveGeodesic(h, X, M);
    }

    double M = 1.989e30; //In solar masses
    double h = 5E-6; //31536/1000

    double previousTime = 0;
    void Update()
    {
        for (int N = 0; N <= 11; N++)
        {
            X = solveGeodesic(h, X, M);

            //UnityEngine.Debug.Log(h);
            //UnityEngine.Debug.Log((X[0] - previousTime)/h);
            h /= (X[0] - previousTime) / h;
            previousTime = X[0];

            //Vector3 v = new Vector3((float)X[3], (float)X[4], (float)X[5]);
            //UnityEngine.Debug.Log(v);
            Tensor XCart = polToCart(X) / scale;
            //Debug.Log(XCart);
            Vector3 finalPos = new Vector3((float)XCart[1], (float)XCart[2], 0);
            //UnityEngine.Debug.Log(finalPos);
            body.localPosition = finalPos;

            //UnityEngine.Debug.Log("solveGeodesic() execution time: " + stopwatch.ElapsedMilliseconds + " ms");
        }
    }
}
