using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

public class ShaderEmitter
{
    //[UnityEditor.Callbacks.DidReloadScripts]
    //[MenuItem("PathTracer/EmitPBRShader %g")]
    public static void SavePBRPathTracer()
    {
        string text = System.IO.File.ReadAllText(Application.dataPath + "/PathTracer_PBR_Template.c"); 
        System.IO.StreamWriter file = new System.IO.StreamWriter(Application.dataPath + "/PathTracer_PBR_Generated.shader");
        file.WriteLine(text);

        file.Close();
        Debug.Log("Shader Generated");
    }




     

}
