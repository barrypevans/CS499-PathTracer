Shader "PathTracer_PBR"
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
			#include "Lighting.cginc"



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
				float4 albedo;
				float metallic;
				float roughness;
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
			uniform samplerCUBE SkyBox;

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
				h.position = 0;
				h.mat.albedo = 0;
				h.mat.metallic = 0;
				h.mat.roughness = 0;
				if (dot(p.normal, -r.direction) < 0)
					return h;

				float t = -dot(r.origin - p.position, p.normal) / dot(p.normal, r.direction);
				if (t < 0)return h;
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
				h.position = 0;
				h.mat.albedo = 0;
				h.mat.metallic = 0;
				h.mat.roughness = 0;
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
				nearest.position = 0;
				nearest.mat.albedo = 0;
				nearest.mat.metallic = 0;
				nearest.mat.roughness = 0;
				float nearestDist = 10000000;

				for (int i = 0; i < SphereBufferLength; i++) {

					Hit h = SphereIntersect(r, SphereBuffer[i]);

					if (h.didhit == 0) continue;

					float dist = distance(h.position, CameraPos);
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



			float3 radiance(float3 position) {
				Ray shadowRay;
				shadowRay.direction = _WorldSpaceLightPos0.xyz;
				shadowRay.origin = position + shadowRay.direction*.0001f;

				Hit h = Intersect(shadowRay);
				if (h.didhit == 0)
					return _LightColor0;
				return float3(0, 0, 0);
			}

			//PBR Functions
			float3 fresnelSchlick(float cosTheta, float3 F0)
			{
				return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
			}

			float DistributionGGX(float3 N, float3 H, float roughness)
			{
				float a = roughness*roughness;
				float a2 = a*a;
				float NdotH = max(dot(N, H), 0.0);
				float NdotH2 = NdotH*NdotH;

				float nom = a2;
				float denom = (NdotH2 * (a2 - 1.0) + 1.0);
				denom = 3.141592f * denom * denom;

				return nom / denom;
			}

			float GeometrySchlickGGX(float NdotV, float roughness)
			{
				float r = (roughness + 1.0);
				float k = (r*r) / 8.0;

				float nom = NdotV;
				float denom = NdotV * (1.0 - k) + k;

				return nom / denom;
			}
			float GeometrySmith(float3 N, float3 V, float3 L, float roughness)
			{
				float NdotV = max(dot(N, V), 0.0);
				float NdotL = max(dot(N, L), 0.0);
				float ggx2 = GeometrySchlickGGX(NdotV, roughness);
				float ggx1 = GeometrySchlickGGX(NdotL, roughness);

				return ggx1 * ggx2;
			}

			float3 fresnelSchlickRoughness(float cosTheta, float3 F0, float roughness)
			{
				return F0 + (max(float3(1.0 - roughness, 1.0 - roughness, 1.0 - roughness), F0) - F0) * pow(1.0 - cosTheta, 5.0);
			}

			//pbr functions end

			//hitpoint and incoming ray
			float3 GetDirectLight(Hit h, Ray ray) {
				float3 N = normalize(h.normal);
				float3 V = normalize(ray.origin - h.position);
				float3 L = normalize(_WorldSpaceLightPos0.xyz);
				float3 H = normalize(V + L);
				float3 rad = radiance(h.position);

				float3 F0 = float3(.04, .04, .04);
				F0 = lerp(F0, h.mat.albedo, h.mat.metallic);
				float3 F = fresnelSchlick(max(dot(H, V), 0.0), F0);
				float NDF = DistributionGGX(N, H, h.mat.roughness);
				float G = GeometrySmith(N, V, L, h.mat.roughness);

				float3 nominator = NDF * G * F;
				float denominator = 4 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.001;
				float3 specular = nominator / denominator;

				float3 kS = F;
				float3 kD = float3(1.0, 1.0, 1.0) - kS;

				kD *= 1.0 - h.mat.metallic;


				float NdotL = max(dot(N, L), 0.0);
				return (kD * h.mat.albedo / 3.141592f + specular) *  rad * NdotL;
			}


			fixed4 frag(v2f i) : SV_Target
			{
				Ray ray = GetRay(CameraPos,i.uv);
			ray.direction = mul((float3x3)LookMatrix, ray.direction);
			Hit h;

			h = Intersect(ray);

			if ((h.didhit) == 0)
				return texCUBE(SkyBox, ray.direction);


			float3 finalColor = float3(0,0,0);
			float3 colorAccum = float3(1,1,1);




				float3 color = float3(0,0,0);

				color = GetDirectLight(h,ray);

				float3 albedo = h.mat.albedo;

				float3 F0 = float3(.04, .04, .04);
				F0 = lerp(F0, h.mat.albedo, h.mat.metallic);
				float3 kS = fresnelSchlickRoughness(max(dot(h.normal, normalize(ray.origin - h.position)), 0.0), F0,h.mat.roughness);
				float3 kD = 1.0 - kS;

				float3 irradiance=float3(.03,.03,.03)*2;


				float3 diff = irradiance *albedo;
				float ambient =kD* diff;


				finalColor += max(0.001f, color)+ ambient;
				 

			finalColor = finalColor / (finalColor + float3(1, 1, 1));
			finalColor = pow(finalColor, float3(1.0 / 2.2, 1.0 / 2.2, 1.0 / 2.2));
			return float4(lerp(finalColor, tex2D(_MainTex,i.uv).xyz, FrameBias), 1);
			//return float4 (finalColor, 1);

			}
			ENDCG
		}
	}
}
