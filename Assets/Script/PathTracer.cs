using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(Camera))]
public class PathTracer : MonoBehaviour {

    private Material _mat;
    private RenderTexture _frameBuffer;
    private RenderTexture _lastFrameBuffer;
    private Camera _cam;

    private ComputeBuffer _sphereBuffer;
    private ComputeBuffer _planeBuffer;
    [SerializeField]
    private Sphere[] _sphereList;
    [SerializeField]
    private Plane[] _planeList;

    public float _frameCount = 0;

    private void Start()
    {
        _cam = GetComponent<Camera>();
        _mat = new Material(Shader.Find("PathTracer"));
        _frameBuffer = new RenderTexture(Screen.width, Screen.height, 16, RenderTextureFormat.ARGBFloat);
        _lastFrameBuffer = new RenderTexture(Screen.width, Screen.height, 16, RenderTextureFormat.ARGBFloat);
        
    }

    private void OnValidate()
    {
        _frameCount = 0;
    }

    private void Update()
    {

        _mat.SetFloat("FOV", _cam.fieldOfView * .5f * Mathf.Deg2Rad);
        if (transform.hasChanged)
        {
            _frameCount = 0;
            transform.hasChanged = false;
        }
    }

    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        if (_sphereList.Length > 0)
        {
            _sphereBuffer = new ComputeBuffer(_sphereList.Length, 36);
            _sphereBuffer.SetData(_sphereList);
            _mat.SetBuffer("SphereBuffer", _sphereBuffer);
            _mat.SetInt("SphereBufferLength", _sphereList.Length);
            //_sphereBuffer.Dispose();
            // _sphereBuffer.Release();
        }

        if (_planeList.Length > 0)
        {
            _planeBuffer = new ComputeBuffer(_planeList.Length, 44);
            _planeBuffer.SetData(_planeList);
            _mat.SetBuffer("PlaneBuffer", _planeBuffer);
            _mat.SetInt("PlaneBufferLength", _planeList.Length);
            // _planeBuffer.Dispose();
            // _planeBuffer.Release();
        }

        _mat.SetFloat("FrameBias", _frameCount / (_frameCount + 1.0f));
        _frameCount++;

        for (int i = 0; i < 6; i++) 
        _mat.SetFloat("rand"+i, Random.Range(0.0f, 1.0f));
        _mat.SetVector("CameraPos", _cam.transform.position);
        _mat.SetMatrix("LookMatrix", new Matrix4x4((Vector4)transform.right, (Vector4)transform.up, (Vector4)transform.forward,Vector4.zero));

        Graphics.Blit(_frameBuffer, _lastFrameBuffer);
        Graphics.Blit(_lastFrameBuffer,_frameBuffer, _mat);
        Graphics.Blit(_frameBuffer,destination);
    }

    [System.Serializable]
    public struct PTMaterial
    {
        public Vector4 color;
        public float surfaceType;
    };

    [System.Serializable]
    public struct Sphere
    {
        public Vector3 position;
        public float radius;
        public PTMaterial mat;
    };
    [System.Serializable]
    public struct Plane
    {
        public Vector3 position;
        public Vector3 normal;
        public PTMaterial mat;
    };


}
