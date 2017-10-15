using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Octree {

    private ComputeShader _computeShader;
    private int _kernel;
    private ComputeBuffer _TestBuffer;

    public int MaxDepth;


    public Octree(int depth)
    {
        MaxDepth = depth;
        _computeShader = Resources.Load<ComputeShader>("Octree");
        _kernel = _computeShader.FindKernel("Add");
        int bucketLength = (int)(Mathf.Pow((float)8, (float)MaxDepth + 1f)-1)/7;
        _TestBuffer = new ComputeBuffer(bucketLength, 4);
        int[] result = new int[bucketLength];
        _TestBuffer.SetData(result);
        _computeShader.SetBuffer(_kernel, "TestBuffer", _TestBuffer);
        _computeShader.Dispatch(_kernel, 4, 1, 1);
        _TestBuffer.GetData(result);
        Debug.Log(result[0]);
    }

    public void Test()
    {
       
      
    }

    public void Insert(Vector3[] Verts, int[] Indices)
    {
        
    }


}

