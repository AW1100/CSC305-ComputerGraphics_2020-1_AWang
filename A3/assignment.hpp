#pragma once

#include "paths.hpp"

#include <exception>
#include <iostream>
#include <string>

#include <atlas/glx/Buffer.hpp>
#include <atlas/glx/Context.hpp>
#include <atlas/glx/ErrorCallback.hpp>
#include <atlas/glx/GLSL.hpp>
#include <atlas/utils/Cameras.hpp>
#include <atlas/utils/LoadObjFile.hpp>
#include <atlas/gui/GUI.hpp>

#include <fmt/printf.h>
#include <magic_enum.hpp>

using namespace atlas;

static constexpr float nearVal{ 1.0f };
static constexpr float farVal{ 10000000000.0f };

static const std::vector<std::string> IncludeDir{ ShaderPath };

struct OpenGLError : std::runtime_error
{
    OpenGLError(const std::string& what_arg) : std::runtime_error(what_arg) {};
    OpenGLError(const char* what_arg) : std::runtime_error(what_arg) {};
};

//class Light
//{
//    void createLight()
//    {
//        GLuint lightVAO;
//        glCreateVertexArrays(1, &lightVAO);
//        glBindVertexArray(lightVAO);
//        glBindBuffer(GL_ARRAY_BUFFER, VBO);
//    };
//
//};

class Triangle
{
public:
    Triangle()
    {
        // allocate the memory to hold the program and shader data
        mProgramHandle = glCreateProgram();
        mVertHandle = glCreateShader(GL_VERTEX_SHADER);
        mFragHandle = glCreateShader(GL_FRAGMENT_SHADER);
        position = 0.0f;
    };

    void loadShaders()
    {
        std::string shaderRoot{ ShaderPath };
        vertexSource =
            glx::readShaderSource(shaderRoot + "triangle.vert", IncludeDir);
        fragmentSource =
            glx::readShaderSource(shaderRoot + "triangle.frag", IncludeDir);

        if (auto result{ glx::compileShader(vertexSource.sourceString, mVertHandle) };
            result)
        {
            throw OpenGLError(*result);
        }

        if (auto result =
            glx::compileShader(fragmentSource.sourceString, mFragHandle);
            result)
        {
            throw OpenGLError(*result);
        }

        // communicate to OpenGL the shaders used to render the Triangle
        glAttachShader(mProgramHandle, mVertHandle);
        glAttachShader(mProgramHandle, mFragHandle);

        if (auto result = glx::linkShaders(mProgramHandle); result)
        {
            throw OpenGLError(*result);
        }

        setupUniformVariables();
    };

    void loadDataToGPU(std::vector<float> const& vertices)
    {
        // create buffer to hold triangle vertex data
        glCreateBuffers(1, &mVbo);
        // allocate and initialize buffer to vertex data
        glNamedBufferStorage(
            mVbo, glx::size<float>(vertices.size()), vertices.data(), 0);

        // create holder for all buffers
        glCreateVertexArrays(1, &mVao);
        // bind vertex buffer to the vertex array
        glVertexArrayVertexBuffer(mVao, 0, mVbo, 0, glx::stride<float>(9));
        

        // enable attributes for the two components of a vertex
        glEnableVertexArrayAttrib(mVao, 0);
        glEnableVertexArrayAttrib(mVao, 1);
        glEnableVertexArrayAttrib(mVao, 2);

        // specify to OpenGL how the vertices and colors are laid out in the buffer
        glVertexArrayAttribFormat(
            mVao, 0, 3, GL_FLOAT, GL_FALSE, glx::relativeOffset<float>(0));
        glVertexArrayAttribFormat(
            mVao, 1, 3, GL_FLOAT, GL_FALSE, glx::relativeOffset<float>(3));
        glVertexArrayAttribFormat(
            mVao, 2, 3, GL_FLOAT, GL_FALSE, glx::relativeOffset<float>(6));

        // associate the vertex attributes (coordinates and color) to the vertex
        // attribute
        glVertexArrayAttribBinding(mVao, 0, 0);
        glVertexArrayAttribBinding(mVao, 1, 0);
        glVertexArrayAttribBinding(mVao, 2, 0);
    };

    void reloadShaders()
    {
        if (glx::shouldShaderBeReloaded(vertexSource))
        {
            glx::reloadShader(
                mProgramHandle, mVertHandle, vertexSource, IncludeDir);
        }

        if (glx::shouldShaderBeReloaded(fragmentSource))
        {
            glx::reloadShader(
                mProgramHandle, mFragHandle, fragmentSource, IncludeDir);
        }
    };

    void render([[maybe_unused]] bool paused,
        [[maybe_unused]] int width,
        [[maybe_unused]] int height)
    {
        reloadShaders();



        if (!paused)
        {
            //position = static_cast<float>(glfwGetTime()) * 64.0f;
            position += 0.50f;
        }

        auto modelMat{ glm::rotate(math::Matrix4{0.1f},
            glm::radians(position),
            math::Vector{0.0f,1.0f,0.0f}) };


        glm::vec3 viewPoint{ 0.0f,0.0f,0.0f };

        auto viewMat{ glm::lookAt(glm::vec3{0.0f,5.0f,20.0f},		//camera looking
                                    viewPoint,		//camera is
                                    glm::vec3{0.0f,1.0f,0.0f}) };

        auto projMat{ glm::perspective(glm::radians(1.0f),
                                static_cast<float>(width) / height,
                                nearVal,
                                farVal) };

        glm::vec3 pointLightPos{ 0.0f,30.0f,-30.0f };
        

        // tell OpenGL which program object to use to render the Triangle
        glUseProgram(mProgramHandle);

        glUniformMatrix4fv(mUniformProjectionLoc, 1, GL_FALSE, glm::value_ptr(projMat));
        glUniformMatrix4fv(mUniformViewLoc, 1, GL_FALSE, glm::value_ptr(viewMat));
        glUniformMatrix4fv(mUniformModelLoc, 1, GL_FALSE, glm::value_ptr(modelMat));

        glUniform3fv(mUniformLightPosLoc, 1, glm::value_ptr(pointLightPos));
        glUniform3fv(mUniformViewPosLoc, 1, glm::value_ptr(viewPoint));

        // tell OpenGL which vertex array object to use to render the Triangle
        glBindVertexArray(mVao);
        // actually render the Triangle
        glDrawArrays(GL_TRIANGLES, 0, 36);

    };

    void freeGPUData()
    {
        // unwind all the allocations made
        glDeleteVertexArrays(1, &mVao);
        glDeleteBuffers(1, &mVbo);
        glDeleteShader(mFragHandle);
        glDeleteShader(mVertHandle);
        glDeleteProgram(mProgramHandle);
    };

private:
    void setupUniformVariables()
    {
        mUniformModelLoc = glGetUniformLocation(mProgramHandle, "model");
        mUniformViewLoc = glGetUniformLocation(mProgramHandle, "view");
        mUniformProjectionLoc = glGetUniformLocation(mProgramHandle, "projection");

        mUniformLightPosLoc = glGetUniformLocation(mProgramHandle, "lightPos");
        mUniformViewPosLoc = glGetUniformLocation(mProgramHandle, "viewPos");
    };

    

    float position;

    // Vertex buffers.
    GLuint mVao;
    GLuint mVbo;

    // Shader data.
    GLuint mVertHandle;
    GLuint mFragHandle;
    GLuint mProgramHandle;
    glx::ShaderFile vertexSource;
    glx::ShaderFile fragmentSource;

    // Uniform variable data.
    GLuint mUniformModelLoc;
    GLuint mUniformViewLoc;
    GLuint mUniformProjectionLoc;

    GLuint mUniformLightPosLoc;
    GLuint mUniformViewPosLoc;
};

class Program
{
public:
    Program(int width, int height, std::string title) :
        settings{}, callbacks{}, paused{}, mWindow{ nullptr }
    {
        settings.size.width = width;
        settings.size.height = height;
        settings.title = title;

        if (!glx::initializeGLFW(errorCallback))
        {
            throw OpenGLError("Failed to initialize GLFW with error callback");
        }

        mWindow = glx::createGLFWWindow(settings);
        if (mWindow == nullptr)
        {
            throw OpenGLError("Failed to create GLFW Window");
        }

        callbacks.keyPressCallback = [&](int key, int, int action, int) { //[=]
            if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE)
            {
                paused = !paused;
            }

        };

        createGLContext();
    };

    void run(Triangle& tri)
    {
        glEnable(GL_DEPTH_TEST);

        ImGui::CreateContext();
        ImGui::StyleColorsDark();

        atlas::gui::GuiWindowData windowData;
        atlas::gui::GuiRenderData renderData;

        gui::initializeGuiWindowData(windowData);
        gui::initializeGuiRenderData(renderData);
        gui::setGuiWindow(windowData, mWindow);

        while (!glfwWindowShouldClose(mWindow))
        {
            int width;
            int height;

            glfwGetFramebufferSize(mWindow, &width, &height);
            // setup the view to be the window's size
            glViewport(0, 0, width, height);
            // tell OpenGL the what color to clear the screen to
            glClearColor(0, 0, 0, 1);
            // actually clear the screen
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            //glClear(GL_COLOR_BUFFER_BIT);
            tri.render(paused, width, height);

            gui::newFrame(windowData);

            ImGui::Begin("Filter controls");

            if (ImGui::Button("Pause"))
            {
                paused = !paused;
            }

            //ImGui::ShowDemoWindow();
            ImGui::End();
            ImGui::Render();
            gui::endFrame(windowData, renderData);

            glfwSwapBuffers(mWindow);
            glfwPollEvents();
        }
    };

    void freeGPUData()
    {
        glx::destroyGLFWWindow(mWindow);
        glx::terminateGLFW();
    };

private:
    static void errorCallback(int code, char const* message)
    {
        fmt::print("error ({}): {}\n", code, message);
    }

    void createGLContext()
    {
        using namespace magic_enum::bitwise_operators;

        glx::bindWindowCallbacks(mWindow, callbacks);
        glfwMakeContextCurrent(mWindow);
        glfwSwapInterval(1);

        if (!glx::createGLContext(mWindow, settings.version))
        {
            throw OpenGLError("Failed to create OpenGL context");
        }

        glx::initializeGLCallback(glx::ErrorSource::All,
            glx::ErrorType::All,
            glx::ErrorSeverity::High |
            glx::ErrorSeverity::Medium);
    };

    GLFWwindow* mWindow;
    glx::WindowSettings settings;
    glx::WindowCallbacks callbacks;

    bool paused;
};
