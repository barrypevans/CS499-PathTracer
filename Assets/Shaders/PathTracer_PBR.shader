Shader "PathTracer_PBR"
{/*
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
				float4 emission;
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

			struct Triangle {
				float3 position;
				float3 normal;
				Material mat;
			};

			struct ILInfo {
				Ray ray;
				Hit hit;
			};


			uniform StructuredBuffer<Triangle> TriangleBuffer;

			uniform int TriangleBufferLength;

			uniform sampler2D _MainTex;
			uniform float4 _MainTex_TexelSize;
			uniform float FOV;
			uniform float FrameBias;
			uniform samplerCUBE SkyBox;

			uniform float3 CameraPos;

			uniform float4x4 LookMatrix;
			uniform uint FrameCount;

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

			Hit TriangleIntersect(Ray r, Triangle p) {

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

				return nearest;
			}



			float3 radiance(float3 position) {
				Ray shadowRay;
				shadowRay.direction = _WorldSpaceLightPos0.xyz;
				shadowRay.origin = position;

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

			float GeometrySchlickGGX_IBL(float NdotV, float roughness)
			{
				float a = roughness;
				float k = (a * a) / 2.0;

				float nom = NdotV;
				float denom = NdotV * (1.0 - k) + k;

				return nom / denom;
			}
			// ----------------------------------------------------------------------------
			float GeometrySmith_IBL(float3 N, float3 V, float3 L, float roughness)
			{
				float NdotV = max(dot(N, V), 0.0);
				float NdotL = max(dot(N, L), 0.0);
				float ggx2 = GeometrySchlickGGX_IBL(NdotV, roughness);
				float ggx1 = GeometrySchlickGGX_IBL(NdotL, roughness);

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
				float3 tangent = (cross(up, N));
				tangent /= length(tangent);
				float3 bitangent = cross(N, tangent);

				float3 sampleVec = tangent * H.x + bitangent * H.y + N * H.z;
				sampleVec /= length(sampleVec);
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
			float2 Hammersley(uint i)
			{
				return float2((rand0* 300.0f) / 300.0f, RadicalInverse_VdC(i));
			}



			//pbr functions end

			float4 GetSky(float3 dir) {
				float4 tex = texCUBE(SkyBox, dir);
				half3 c = DecodeHDR(tex, unity_SpecCube0_HDR);
				return tex;
			}

	 

			 
			float3 GetDirectLight(Hit h, Ray ray) {
				if (h.didhit == 0) {
					return GetSky(ray.direction);
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
				float G = GeometrySmith_IBL(N, V, L, h.mat.roughness);

				float3 nominator = NDF * G * F;
				float denominator = 4 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.001;
				float3 specular = nominator / denominator;

				float3 kS = F;
				float3 kD = float3(1.0, 1.0, 1.0) - kS;

				kD *= 1.0 - h.mat.metallic;


				float NdotL = max(dot(N, L), 0.0);

				float3 lo = (kD * h.mat.albedo / 3.141592f + specular) *  rad * NdotL;
					
				//emission next event estimation
				for (int i = 0; i < SphereBufferLength; i++) {
						if (length(SphereBuffer[i].mat.emission) > 0)continue;

						Ray shadowRay;

						float3 rpoint = (float3(rand1, rand3, rand0) - .5) * 2.0f;

						rpoint /= length(rpoint);
						rpoint *= SphereBuffer[i].radius;
						rpoint += SphereBuffer[i].position;

						shadowRay.direction = normalize(rpoint - h.position);
						shadowRay.origin = h.position;

						Hit shadowh = Intersect(shadowRay);
						if (shadowh.didhit == 0)continue;
						if (length(shadowh.mat.emission) > 0)continue;

						L = normalize(shadowRay.direction);
						H = normalize(V + L);

						float dist = length(rpoint - h.position);
						float attenuation = 1.0 / (dist * dist);
						 rad = shadowh.mat.emission*attenuation;

						 F0 = float3(.04, .04, .04);
						F0 = lerp(F0, h.mat.albedo, h.mat.metallic);
						 F = fresnelSchlick(max(dot(H, V), 0.0), F0);
						 NDF = DistributionGGX(N, H, h.mat.roughness);
						 G = GeometrySmith_IBL(N, V, L, h.mat.roughness);

						 nominator = NDF * G * F;
						 denominator = 4 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.001;
						 specular = nominator / denominator;

						 kS = F;
						 kD = float3(1.0, 1.0, 1.0) - kS;

						kD *= 1.0 - h.mat.metallic;


						 NdotL = max(dot(N, L), 0.0);


						lo += (kD * h.mat.albedo / 3.141592f + specular) *  rad * NdotL;
					}


				return lo;
			}


			fixed4 frag(v2f i) : SV_Target
			{
				Ray ray = GetRay(CameraPos,i.uv);
			ray.direction = mul((float3x3)LookMatrix, ray.direction);
			Hit h;

			h = Intersect(ray);

			if ((h.didhit) == 0)
				return GetSky(ray.direction);


			float3 finalColor = float3(0,0,0);

			ILInfo Stack[STACKSIZE];
			float3 FinalStack[STACKSIZE];

			for (int l = 0; l < STACKSIZE; l++) FinalStack[l] = float3(0, 0, 0);


			Stack[0].ray = ray;
			Stack[0].hit = h;

			for (int k = 1; k < STACKSIZE; k++) {



				int parentIndex = ceil((float(k) - 1) / 2.0f);
				if (k % 2 == 0) {
					Ray newRay;
					newRay.direction = sampleUnitHemisphereNormal(Stack[parentIndex].hit.normal);
					newRay.origin = Stack[parentIndex].hit.position;

					h = Intersect(newRay);
					Stack[k].hit = h;
					Stack[k].ray = newRay;
					if (BOUNCES == ceil(log(k + 2) / log(2))) {
						FinalStack[k] = GetDirectLight(Stack[k].hit, Stack[k].ray) + Stack[k].hit.mat.emission;
						continue;
					}
				}
				else {
					Ray newRay;
					float3 N = normalize(Stack[parentIndex].hit.normal);
					float3 V = normalize(Stack[parentIndex].ray.origin - Stack[parentIndex].hit.position);

					float2 Xi = Hammersley(uint(FrameCount));
					float3 H = ImportanceSampleGGX(Xi, N, Stack[parentIndex].hit.mat.roughness);
					float3 L = normalize(-reflect(V,H));
					newRay.direction = L;
					newRay.origin = Stack[parentIndex].hit.position;

					h = Intersect(newRay);
					Stack[k].hit = h;
					Stack[k].ray = newRay;
					if (BOUNCES == ceil(log(k + 2) / log(2))) {
						FinalStack[k] = GetDirectLight(Stack[k].hit, Stack[k].ray) + Stack[k].hit.mat.emission;
						continue;
					}
				}
			}


			for (int j = STACKSIZE - pow(2,BOUNCES) - 1; j >= 0; j--) {
				if (Stack[j].hit.didhit == 0 && j != 0) {
					if (Stack[ceil((float(j) - 1) / 2.0f)].hit.didhit == 0) {
						FinalStack[j] = float3(0, 0, 0);
						continue;
					}
				}
				float3 N = normalize(Stack[j].hit.normal);
				float3 V = normalize(Stack[j].ray.origin - Stack[j].hit.position);
				float3 l = FinalStack[(j * 2) + 1];
				float3 r = FinalStack[(j * 2) + 2];
				float3 irradiance = r;
				irradiance /= 2;
				float3 F0 = float3(.04, .04, .04);
				F0 = lerp(F0, Stack[j].hit.mat.albedo, Stack[j].hit.mat.metallic);
				float3 kS = fresnelSchlickRoughness(max(dot(N, V), 0.0), F0, Stack[j].hit.mat.roughness);
				float3 kD = 1.0 - kS;
				float3 diffuse = irradiance * Stack[j].hit.mat.albedo;
				float3 ambient = (kD * diffuse + kS*l);
				if (j % 2 == 0)
					FinalStack[j] = ambient + GetDirectLight(Stack[j].hit, Stack[j].ray) + Stack[j].hit.mat.emission;
				else
					FinalStack[j] = ambient + GetDirectLight(Stack[j].hit, Stack[j].ray) + Stack[j].hit.mat.emission;

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
	}*/
}
