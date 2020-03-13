#version 450 core

in vec3 vertexColour;
in vec3 Normal;
in vec3 FragPos;

out vec4 fragColour;

uniform vec3 lightPos;
uniform vec3 viewPos;

void main()
{
    float ambientStrength = 0.1;
    vec3 ambient = ambientStrength * vec3(1.0,1.0,1.0);

    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * vec3(1.0,1.0,1.0);

    float specularStrength = 1.5;
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);  
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 64);
    vec3 specular = specularStrength * spec * vec3(1.0,1.0,1.0);  

    vec3 result = (ambient + diffuse + specular) * vertexColour;
    fragColour = vec4(result, 1.0);
}
