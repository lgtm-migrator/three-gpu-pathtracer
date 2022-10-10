export const newSampling = /* glsl */`
	float GTR2Aniso(float NDotH, float HDotX, float HDotY, float ax, float ay)
	{
		float a = HDotX / ax;
		float b = HDotY / ay;
		float c = a * a + b * b + NDotH * NDotH;
		return 1.0 / (PI * ax * ay * c * c);
	}

	vec3 SampleGGXVNDF(vec3 V, float ax, float ay, float r1, float r2)
	{
		vec3 Vh = normalize(vec3(ax * V.x, ay * V.y, V.z));

		float lensq = Vh.x * Vh.x + Vh.y * Vh.y;
		vec3 T1 = lensq > 0.0 ? vec3(-Vh.y, Vh.x, 0) * inversesqrt(lensq) : vec3(1, 0, 0);
		vec3 T2 = cross(Vh, T1);

		float r = sqrt(r1);
		float phi = 2.0 * PI * r2;
		float t1 = r * cos(phi);
		float t2 = r * sin(phi);
		float s = 0.5 * (1.0 + Vh.z);
		t2 = (1.0 - s) * sqrt(1.0 - t1 * t1) + s * t2;

		vec3 Nh = t1 * T1 + t2 * T2 + sqrt(max(0.0, 1.0 - t1 * t1 - t2 * t2)) * Vh;

		return normalize(vec3(ax * Nh.x, ay * Nh.y, max(0.0, Nh.z)));
	}

	float SchlickFresnel(float u)
	{
		float m = clamp(1.0 - u, 0.0, 1.0);
		float m2 = m * m;
		return m2 * m2 * m;
	}

	float DielectricFresnel(float cosThetaI, float eta)
	{
		float sinThetaTSq = eta * eta * (1.0f - cosThetaI * cosThetaI);

		// Total internal reflection
		if (sinThetaTSq > 1.0)
			return 1.0;

		float cosThetaT = sqrt(max(1.0 - sinThetaTSq, 0.0));

		float rs = (eta * cosThetaT - cosThetaI) / (eta * cosThetaT + cosThetaI);
		float rp = (eta * cosThetaI - cosThetaT) / (eta * cosThetaI + cosThetaT);

		return 0.5f * (rs * rs + rp * rp);
	}

	float Luminance(vec3 c)
	{
		return 0.212671 * c.x + 0.715160 * c.y + 0.072169 * c.z;
	}


	float DisneyFresnel(float metalness, float eta, float LDotH, float VDotH) {

		float metallicFresnel = SchlickFresnel(LDotH);
		float dielectricFresnel = DielectricFresnel(abs(VDotH), eta);
		return mix(dielectricFresnel, metallicFresnel, metalness);

	}

	float SmithGAniso(float NDotV, float VDotX, float VDotY, float ax, float ay)
	{
		float a = VDotX * ax;
		float b = VDotY * ay;
		float c = NDotV;
		return (2.0 * NDotV) / (NDotV + sqrt(a * a + b * b + c * c));
	}

	void GetSpecColor(SurfaceRec surf, float eta, out vec3 specCol, out vec3 sheenCol)
	{
		float lum = Luminance(surf.color);
		vec3 ctint = lum > 0.0 ? surf.color / lum : vec3(1.0f);
		float F0 = (1.0 - eta) / (1.0 + eta);
		specCol = mix(F0 * F0 * mix(vec3(1.0), ctint, surf.specularColor), surf.color, surf.metalness);
		sheenCol = mix(vec3(1.0), ctint, surf.sheenColor);
	}


`;
