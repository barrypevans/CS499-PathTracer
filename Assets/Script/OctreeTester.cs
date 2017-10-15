using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class OctreeTester : MonoBehaviour {

    Octree octree;

	// Use this for initialization
	void Start () {
        octree = new Octree(3);
        octree.Test();
    }
	
	// Update is called once per frame
	void Update () {
		
	}
}
