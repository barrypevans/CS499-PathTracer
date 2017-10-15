Shader "PathTracer_PBR_diffuseindirect"
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
			#define	BOUNCES 3
			#define	STACKSIZE pow(2,BOUNCES+1)-1

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

			struct ILInfo {
				Ray ray;
				Hit hit;
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

			float3 ImportanceSampleGGX(float2 Xi, float3 N, float roughness)
			{
				float a = roughness*roughness;

				float phi = 2.0 * 3.14159 * Xi.x;
				float cosTheta = sqrt((1.0 - Xi.y) / (1.0 + (a*a - 1.0) * Xi.y));
				float sinTheta = sqrt(1.0 - cosTheta*cosTheta);

				// from spherical coordinates to cartesian coordinates
				float3 H;
				H.x = cos(phi) * sinTheta;
				H.y = sin(phi) * sinTheta;
				H.z = cosTheta;

				// from tangent-space vector to world-space sample vector
				float3 up = abs(N.z) < 0.999 ? float3(0.0, 0.0, 1.0) : float3(1.0, 0.0, 0.0);
				float3 tangent = normalize(cross(up, N));
				float3 bitangent = cross(N, tangent);

				float3 sampleVec = tangent * H.x + bitangent * H.y + N * H.z;
				return normalize(sampleVec);
			}

			float RadicalInverse_VdC(uint bits)
			{
				bits = (bits << 16u) | (bits >> 16u);
				bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
				bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
				bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
				bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
				return float(bits) * 2.3283064365386963e-10; // / 0x100000000
			}
			// ----------------------------------------------------------------------------
			float2 Hammersley(uint i, uint N)
			{
				return float2(float(i) / float(N), RadicalInverse_VdC(i));
			}

			//pbr functions end

			//hitpoint and incoming ray
			float3 GetDirectLight(Hit h, Ray ray) {
				if (h.didhit == 0) {
					return texCUBE(SkyBox, ray.direction);	
				}
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

			ILInfo Stack[STACKSIZE];
			float3 FinalStack[STACKSIZE];

			for (int l = 0; l < STACKSIZE; l++) FinalStack[l] = float3(0, 0, 0);
			

			Stack[0].ray = ray;
			Stack[0].hit = h;

			for (int k = 1; k < STACKSIZE; k++) {

				int parentIndex = ceil((float(k) - 1) / 2.0f);
			
					Ray newRay;
					newRay.direction = sampleUnitHemisphereNormal(Stack[parentIndex].hit.normal);
					newRay.origin = Stack[parentIndex].hit.position;

					h = Intersect(newRay);
					Stack[k].hit = h;
					Stack[k].ray = newRay;

				if (BOUNCES == ceil(log(k + 2) / log(2))) {
					FinalStack[k] = GetDirectLight(Stack[k].hit, Stack[k].ray);
				}

			}


			for (int j = STACKSIZE - pow(2,BOUNCES) - 1; j >= 0; j--) {
				float3 N = normalize(Stack[j].hit.normal);
				float3 V = normalize(Stack[j].ray.origin - Stack[j].hit.position);
				float3 l = FinalStack[(j * 2) + 1];
				float3 r = FinalStack[(j * 2) + 2];
				float3 irradiance = l + r;
				irradiance /= 2;
				float3 F0 = float3(.04, .04, .04);
				F0 = lerp(F0, Stack[j].hit.mat.albedo, Stack[j].hit.mat.metallic);
				float3 kS = fresnelSchlickRoughness(max(dot(N, V), 0.0), F0, Stack[j].hit.mat.roughness);
				float3 kD = 1.0 - kS;
				float3 diffuse = irradiance * Stack[j].hit.mat.albedo;
				float3 ambient = (kD * diffuse) ;
				FinalStack[j] = ambient + GetDirectLight(Stack[j].hit, Stack[j].ray);
				
			}

			float3 color = float3(0, 0, 0);

			color = GetDirectLight(h, ray);

			float3 albedo = h.mat.albedo;



			finalColor += max(0.001f, FinalStack[0]);


			finalColor = finalColor / (finalColor + float3(1, 1, 1));
			finalColor = pow(finalColor, float3(1.0 / 2.2, 1.0 / 2.2, 1.0 / 2.2));
			return float4(lerp(finalColor, tex2D(_MainTex,i.uv).xyz, FrameBias), 1);
			//return float4 (finalColor, 1);

			}
			ENDCG
		}
	}
}
