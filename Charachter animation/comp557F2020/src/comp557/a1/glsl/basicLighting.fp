#version 400

//**** Marwen Kraiem 260955202 ****//

uniform vec3 kd;
uniform vec3 ks;
uniform vec3 ka;
   
uniform vec3 lightDir_1;
uniform vec3 lightColor_1;

uniform vec3 lightDir_2;
uniform vec3 lightColor_2;

uniform vec3 lightDir_3;
uniform vec3 lightColor_3;

uniform float Shininess;
uniform vec3 ambientLightColor;

in vec3 normalForFP;
in vec3 FragPos;

out vec4 fragColor;

// TODO: Objective 7, GLSL lighting

void main(void) {

	vec3 LambertianLight;//diffuse
	vec3 ambientLight  ;//ambient
	vec3 reflectedLight;//specular
	vec3 scatteredLight;
	vec3 v = normalize(vec3(0, 0, 2.5f)-FragPos);

	float specular;
	float diffuse;

	//ambientLight	
	ambientLight =  ka * ambientLightColor;

	//***** first direction *****//	
	diffuse = max(0.0f, dot(normalForFP,lightDir_1));
	specular = max(0.0f, dot(normalForFP,normalize(v+lightDir_1)));
		
	//reflectedLight
	if ( diffuse== 0.0f)
		specular = 0.0f;
	else
		specular = pow(specular, Shininess); 
		

	//LambertianLight
	LambertianLight =  lightColor_1 * diffuse;
	
	//reflectedLight
	reflectedLight = lightColor_1 * specular ;
	
		
	//***** second direction *****//	
	diffuse = max(0.0f, dot(normalForFP,lightDir_2));
	specular = max(0.0f, dot(normalForFP,normalize(v+lightDir_2)));
		
	//reflectedLight
	if ( diffuse== 0.0f)
		specular = 0.0f;
	else
		specular = pow(specular, Shininess); 
		

	//LambertianLight
	LambertianLight +=  lightColor_2*diffuse;
	
	//reflectedLight
	reflectedLight += lightColor_2 * specular ;
	
	
	//***** third direction *****//	
	diffuse = max(0.0f, dot(normalForFP,lightDir_3));
	specular = max(0.0f, dot(normalForFP,normalize(v+lightDir_3)));
		
	//reflectedLight
	if ( diffuse== 0.0f)
		specular = 0.0f;
	else
		specular = pow(specular, Shininess); 
		

	//LambertianLight
	LambertianLight +=  lightColor_3*diffuse;
	//reflectedLight
	reflectedLight += lightColor_3* specular ;	
	
	
	scatteredLight = ambientLight+kd * LambertianLight;
	
	fragColor = vec4( min(scatteredLight+ ks*reflectedLight,vec3(1.0)), 1 );
}