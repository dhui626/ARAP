#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#define FREEGLUT_STATIC
#include "gl_core_3_3.h"
#include <GL/glut.h>
#include <GL/freeglut_ext.h>

#define TW_STATIC
#include <AntTweakBar.h>


#include <ctime>
#include <memory>
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>

#include "objloader.h"
#include "glprogram.h"
#include "MyImage.h"
#include "VAOImage.h"
#include "VAOMesh.h"
#include "trackball.h"

#include "matlab_utils.h"


GLProgram MyMesh::prog, MyMesh::pickProg, MyMesh::pointSetProg;
GLTexture MyMesh::colormapTex;

MyMesh M;
int viewport[4] = { 0, 0, 1280, 960 };
int actPrimType = MyMesh::PE_VERTEX;

int Iteration = 10;

bool showATB = true;
bool use_uniform_weight = true;
bool use_cot_weight = false;

using MatX3f = Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>;
using MatX3i = Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>;
MatX3f meshX;
MatX3i meshF;

void deform_preprocess()
{
    using namespace Eigen;

    //使用初始计算的ARAP参数化结果作为2D网格变形初值

    matlabEval("clear;");

    eigen2matlab("x", meshX.cast<double>());    //顶点
    eigen2matlab("t", meshF.cast<double>());    //邻接关系
    matlabEval("t = double(t+1);");   // change index to matlab 1-based
    
    matlabEval("ARAP_pre;");

    // 返回初始ARAP参数化结果

    Matrix<float, Dynamic, Dynamic, RowMajor> y;
    matlab2eigen("single([uv zeros(v_count,1)])", y, true);
    if (y.cols() > 3)  y = y.leftCols(3);
    if (y.rows() == 0 || y.cols() != 3) return;

    M.upload(y.data(), y.rows(), nullptr, 0, nullptr);

}


void meshDeform()
{
    using namespace Eigen;

    std::vector<int> P2PVtxIds = M.getConstrainVertexIds();
    std::vector<float> P2PDsts = M.getConstrainVertexCoords();

    vector2matlab("P2PVtxIds", P2PVtxIds);
    eigen2matlab("P2PDsts", Map<MatX3f>(P2PDsts.data(), P2PVtxIds.size(), 3));
    matlabEval("P2PVtxIds = double(P2PVtxIds+1);");   // change index to matlab 1-based

    std::vector<int> itr = { Iteration };
    vector2matlab("iteration", itr);
    matlabEval("ARAP_cal;");

    Matrix<float, Dynamic, Dynamic, RowMajor> y;
    matlab2eigen("single([uv zeros(v_count,1)])", y, true);
    if (y.cols() > 3)  y = y.leftCols(3);
    if (y.rows() == 0 || y.cols() != 3) return;

    M.upload(y.data(), y.rows(), nullptr, 0, nullptr);
}

int mousePressButton;
int mouseButtonDown;
int mousePos[2];

bool msaa = true;


void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    if (msaa) glEnable(GL_MULTISAMPLE);
    else glDisable(GL_MULTISAMPLE);

    glViewport(0, 0, viewport[2], viewport[3]);
    M.draw(viewport);

    if (showATB) TwDraw();
    glutSwapBuffers();

    //glFinish();
}

void onKeyboard(unsigned char code, int x, int y)
{
    if (!TwEventKeyboardGLUT(code, x, y)) {
        switch (code) {
        case 17:
            exit(0);
        case 'f':
            glutFullScreenToggle();
            break;
        case ' ':
            showATB = !showATB;
            break;
        }
    }

    glutPostRedisplay();
}

void onMouseButton(int button, int updown, int x, int y)
{
    if (!showATB || !TwEventMouseButtonGLUT(button, updown, x, y)) {
        mousePressButton = button;
        mouseButtonDown = updown;

        if (updown == GLUT_DOWN) {
            if (button == GLUT_LEFT_BUTTON) {
                if (glutGetModifiers()&GLUT_ACTIVE_CTRL) {
                }
                else {
                    int r = M.pick(x, y, viewport, M.PE_VERTEX, M.PO_ADD);
                }
            }
            else if (button == GLUT_RIGHT_BUTTON) {
                M.pick(x, y, viewport, M.PE_VERTEX, M.PO_REMOVE);
            }
        }
        else { // updown == GLUT_UP
            if (button == GLUT_LEFT_BUTTON);
        }

        mousePos[0] = x;
        mousePos[1] = y;
    }
    glutPostRedisplay();
}


void onMouseMove(int x, int y)
{
    if (!showATB || !TwEventMouseMotionGLUT(x, y)) {
        if (mouseButtonDown == GLUT_DOWN) {
            if (mousePressButton == GLUT_MIDDLE_BUTTON) {
                M.moveInScreen(mousePos[0], mousePos[1], x, y, viewport);
            }
            else if (mousePressButton == GLUT_LEFT_BUTTON) {
                if (!M.moveCurrentVertex(x, y, viewport)) {
                    meshDeform();
                }
                else {
                    M.moveInScreen(mousePos[0], mousePos[1], x, y, viewport);
                }
            }
        }
    }

    mousePos[0] = x; mousePos[1] = y;

    glutPostRedisplay();
}


void onMouseWheel(int wheel_number, int direction, int x, int y)
{
    M.mMeshScale *= direction > 0 ? 1.1f : 0.9f;
    glutPostRedisplay();
}

int initGL(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_MULTISAMPLE);
    glutInitWindowSize(960, 960);
    glutInitWindowPosition(200, 50);
    glutCreateWindow(argv[0]);

    // !Load the OpenGL functions. after the opengl context has been created
    if (ogl_LoadFunctions() == ogl_LOAD_FAILED)
        return -1;

    glClearColor(1.f, 1.f, 1.f, 0.f);

    glutReshapeFunc([](int w, int h) { viewport[2] = w; viewport[3] = h; TwWindowSize(w, h); });
    glutDisplayFunc(display);
    glutKeyboardFunc(onKeyboard);
    glutMouseFunc(onMouseButton);
    glutMotionFunc(onMouseMove);
    glutMouseWheelFunc(onMouseWheel);
    glutCloseFunc([]() {exit(0); });
    return 0;
}


void createTweakbar()
{
    TwBar *bar = TwGetBarByName("MeshViewer");
    if (bar)    TwDeleteBar(bar);

    //Create a tweak bar
    bar = TwNewBar("MeshViewer");
    TwDefine(" MeshViewer size='220 150' color='114 83 52' text=dark alpha=128 position='5 5'"); // change default tweak bar size and color

    TwAddVarRO(bar, "#Vertex", TW_TYPE_INT32, &M.nVertex, " group='Mesh View'");
    TwAddVarRO(bar, "#Face", TW_TYPE_INT32, &M.nFace, " group='Mesh View'");

    TwAddVarRW(bar, "Point Size", TW_TYPE_FLOAT, &M.pointSize, " group='Mesh View' ");
    TwAddVarRW(bar, "Vertex Color", TW_TYPE_COLOR4F, M.vertexColor.data(), " group='Mesh View' help='mesh vertex color' ");


    TwAddVarRW(bar, "Edge Width", TW_TYPE_FLOAT, &M.edgeWidth, " group='Mesh View' ");
    TwAddVarRW(bar, "Edge Color", TW_TYPE_COLOR4F, M.edgeColor.data(), " group='Mesh View' help='mesh edge color' ");

    TwAddVarRW(bar, "Face Color", TW_TYPE_COLOR4F, M.faceColor.data(), " group='Mesh View' help='mesh face color' ");

    TwDefine(" MeshViewer/'Mesh View' opened=false ");

    TwAddVarRW(bar, "Iteration", TW_TYPE_INT32, &Iteration, "");
}

int main(int argc, char *argv[])
{
    if (initGL(argc, argv)) {
        fprintf(stderr, "!Failed to initialize OpenGL!Exit...");
        exit(-1);
    }

    MyMesh::buildShaders();

    std::vector<float> x;
    std::vector<int> f;

    const char* meshpath = argc > 1 ? argv[1] : "../obj/Cow_dABF.obj";
    readObj(meshpath, x, f);

    meshX = Eigen::Map<MatX3f>(x.data(), x.size() / 3, 3);
    meshF = Eigen::Map<MatX3i>(f.data(), f.size() / 3, 3);

    M.upload(x.data(), x.size() / 3, f.data(), f.size() / 3, nullptr);

    //////////////////////////////////////////////////////////////////////////
    TwInit(TW_OPENGL_CORE, NULL);
    //Send 'glutGetModifers' function pointer to AntTweakBar;
    //required because the GLUT key event functions do not report key modifiers states.
    TwGLUTModifiersFunc(glutGetModifiers);
    glutSpecialFunc([](int key, int x, int y) { TwEventSpecialGLUT(key, x, y); glutPostRedisplay(); }); // important for special keys like UP/DOWN/LEFT/RIGHT ...
    TwCopyStdStringToClientFunc([](std::string& dst, const std::string& src) {dst = src; });

    createTweakbar();

    //////////////////////////////////////////////////////////////////////////
    atexit([] { TwDeleteAllBars();  TwTerminate(); }); 

    glutTimerFunc(1, [](int) {
        getMatEngine().connect(""); 
        deform_preprocess();
    }, 
        0);


    //////////////////////////////////////////////////////////////////////////
    glutMainLoop();

    return 0;
}
