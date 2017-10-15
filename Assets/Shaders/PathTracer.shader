Shader "PathTracer"
{
	Properties
	{
		_MainTex("Texture", 2D) = "white" {}
	}
		SubShader
	{
		// No culling or depth
		Cull Off ZWrite Off ZTest Always

		Pass
		{

			Tags{ "LightMode" = "ForwardBase" }

			CGPROGRAM
			#pragma vertex vert
			#pragma fragment frag

			#include "UnityCG.cginc"

			struct appdata
			{
				float4 vertex : POSITION;
				float2 uv : TEXCOORD0;
			};

			struct v2f
			{
				float2 uv : TEXCOORD0;
				float4 vertex : SV_POSITION;
			};

			v2f vert(appdata v)
			{
				v2f o;
				o.vertex = UnityObjectToClipPos(v.vertex);
				o.uv = v.uv;
				return o;
			}
			uniform float rand0;
			uniform float rand1;
			uniform float rand2;
			uniform float rand3;
			uniform float rand4;
			uniform float rand5;
	


			float rand(float3 co)
			{
				return frac(sin(dot(co.xyz, float3(12.9898, 78.233, 45.5432))) * 43758.5453);
			}


			float3 sampleUnitHemisphereNormal(float3 n)
			{
				float3 rvec = (normalize(float3(rand(rand0), rand(rand1), rand(rand2))) - .5) * 2;
				float3 rvec2 = (normalize(float3(rand(rand3), rand(rand4), rand(rand5))) - .5) * 2;

				float3 tangent = normalize(rvec - n * dot(rvec, n));
				float3 bitangent = cross(n, tangent);
				float3x3 tbn = float3x3(tangent, bitangent, n);

				return mul(tbn,rvec2);
			}

			struct Material {
				float4 color;
				float surfaceType;
			};

			struct Ray {
				fixed3 direction;
				float3 origin;
			};

			struct Hit {
				fixed3 position;
				float3 normal;
				Material mat;
				float didhit;
			};

			struct Sphere {
				float3 position;
				float radius;
				Material mat;
			};

			struct Plane {
				float3 position;
				float3 normal;
				Material mat;
			};

			uniform StructuredBuffer<Sphere> SphereBuffer;
			uniform StructuredBuffer<Plane> PlaneBuffer;

			uniform int SphereBufferLength;
			uniform int PlaneBufferLength;

			uniform sampler2D _MainTex;
			uniform float4 _MainTex_TexelSize;
			uniform float FOV;
			uniform float FrameBias;

			uniform float3 CameraPos;

			uniform float4x4 LookMatrix;

			Ray GetRay(float3 origin ,float2 uv) {
				uv = (uv - .5) * 2;
				uv.x *= tan(FOV);
				uv.y *= (_MainTex_TexelSize.w / _MainTex_TexelSize.z)*tan(FOV);
				Ray ray;
				ray.direction = fixed3(uv, 1);
				ray.direction /= length(ray.direction);
				ray.origin = origin;
				return ray;
			}

			Hit PlaneIntersect(Ray r, Plane p) {
				Hit h;
				h.didhit = 0;
				h.normal = float3(0, 0, 0);
				if (dot(p.normal, -r.direction) < 0)
					return h;
				
				float t = -dot(r.origin- p.position, p.normal) / dot(p.normal, r.direction);
				//if (t<0)return h;
				h.didhit = 1;
				h.normal = p.normal;
				h.position = r.origin + r.direction*t;
				h.mat = p.mat;
				return h;
			}

			Hit SphereIntersect(Ray r, Sphere s) {
				Hit h;
				h.didhit = 0;
				h.normal = float3(0, 0, 0);
				float A = dot(r.direction, r.direction);
				float3 OtoP = (r.origin - s.position);
				float B = 2.0f * dot(OtoP, r.direction);
				float C = dot(OtoP, OtoP) - (s.radius*s.radius);
				float disc = (B*B - 4 * A * C);
				if (disc < 0) return h;
				if (dot((s.position - r.origin), r.direction) < 0) return h;
				h.didhit = 1;
				disc = sqrt(disc);
				float t1 = (-B + disc) / (2 * A);
				float t2 = (-B - disc) / (2 * A);
				float t = t1 < t2 ? t1 : t2;
				h.mat = s.mat;
				h.position = r.origin + r.direction * t;
				h.normal = (h.position - s.position) / length(h.position - s.position);
				return h;
			}


			Hit Intersect(Ray r) {

				Hit nearest;
				nearest.didhit = 0;
				nearest.normal = float3(0, 0, 0);
				float nearestDist = 10000000;

				for (int i = 0; i < SphereBufferLength; i++) {

					Hit h = SphereIntersect(r, SphereBuffer[i]);

					if (h.didhit == 0) continue;

					float dist = distance(h.position, float3(0, 0, 0));
					if (dist < nearestDist) {
						nearest = h;
						nearestDist = dist;
					}
				}

				for (int i = 0; i < PlaneBufferLength; i++) {

					Hit h = PlaneIntersect(r, PlaneBuffer[i]);

					if (h.didhit == 0) continue;

					float dist = distance(h.position, float3(0, 0, 0));
					if (dist < nearestDist) {
						nearest = h;
						nearestDist = dist;
					}
				}
				return nearest;
			}



			float3 radiance(float3 position, float3 normal) {
				Ray shadowRay;
				shadowRay.direction = _WorldSpaceLightPos0.xyz;
				shadowRay.origin = position + shadowRay.direction*.0001f;
				
				Hit h = Intersect(shadowRay);
				if (h.didhit == 0)
					return dot(normal, _WorldSpaceLightPos0.xyz);
				return float3(0, 0, 0);
			}



			fixed4 frag(v2f i) : SV_Target
			{
				Ray ray = GetRay(CameraPos,i.uv + (float2(rand(rand0), rand(rand1)) - .5) * 2 * _MainTex_TexelSize.xy);
			ray.direction = mul((float3x3)LookMatrix, ray.direction);
			Hit h = Intersect(ray);

			if ((h.didhit) == 0)
				return float4(0, 0, 0, 1);


			float3 finalColor = float3(0,0,0);
			float3 colorAccum = float3(1,1,1);

			for (int j = 0; j < 15; j++) {
				h = Intersect(ray);

				float3 color = h.mat.color;
				

				
				ray.origin = h.position;
				if (h.mat.surfaceType == 1) {
					float3 light = pow(radiance(h.position, reflect(ray.direction, h.normal)),5);
					ray.direction = reflect(ray.direction, h.normal);
					colorAccum *= color;
					finalColor += colorAccum*light;
				}
				else {
					float3 light = radiance(h.position, h.normal);
					ray.direction = sampleUnitHemisphereNormal(h.normal);
					colorAccum *= color;
					finalColor += max(0,colorAccum)*light;
				}
			}

			return float4(lerp(finalColor, tex2D(_MainTex,i.uv).xyz, FrameBias), 1);
			//return float4 (finalColor, 1);

			}
			ENDCG
		}
	}
}
