#include <iostream>
#include <opencv2/opencv.hpp>
#include <utility>

#include "global.hpp"
#include "rasterizer.hpp"
#include "Triangle.hpp"
#include "Shader.hpp"
#include "Texture.hpp"
#include "OBJ_Loader.h"

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos) {
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0],
            0, 1, 0, -eye_pos[1],
            0, 0, 1, -eye_pos[2],
            0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float angle) {
    Eigen::Matrix4f rotation;
    angle = angle * MY_PI / 180.f;
    rotation << cos(angle), 0, sin(angle), 0,
            0, 1, 0, 0,
            -sin(angle), 0, cos(angle), 0,
            0, 0, 0, 1;

    Eigen::Matrix4f scale;
    scale << 2.5, 0, 0, 0,
            0, 2.5, 0, 0,
            0, 0, 2.5, 0,
            0, 0, 0, 1;

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

    return translate * rotation * scale;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar) {
    Eigen::Matrix4f projection;

    projection << 1 / (aspect_ratio * tan(eye_fov / 2)), 0, 0, 0,
            0, 1 / tan(eye_fov / 2), 0, 0,
            0, 0, -(zNear + zFar) / (zNear - zFar), -2 * zNear * zFar / (zNear - zFar),
            0, 0, -1, 0;

    return projection;
}

Eigen::Vector3f vertex_shader(const vertex_shader_payload &payload) {
    return payload.position;
}

Eigen::Vector3f normal_fragment_shader(const fragment_shader_payload &payload) {
    Eigen::Vector3f return_color = (payload.normal.head<3>().normalized() + Eigen::Vector3f(1.0f, 1.0f, 1.0f)) / 2.f;
    Eigen::Vector3f result;
    result << return_color.x() * 255, return_color.y() * 255, return_color.z() * 255;
    return result;
}

static Eigen::Vector3f reflect(const Eigen::Vector3f &vec, const Eigen::Vector3f &axis) {
    auto costheta = vec.dot(axis);
    return (2 * costheta * axis - vec).normalized();
}

struct light {
    Eigen::Vector3f position;
    Eigen::Vector3f intensity;
};

Eigen::Vector3f texture_fragment_shader(const fragment_shader_payload &payload) {
    Eigen::Vector3f return_color = {0, 0, 0};
    if (payload.texture) {
        return_color = payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y());
    }
    Eigen::Vector3f texture_color;
    texture_color << return_color.x(), return_color.y(), return_color.z();

    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = texture_color / 255.f;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20,  20,  20},
                    {500, 500, 500}};
    auto l2 = light{{-20, 20,  0},
                    {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};
    Eigen::Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Eigen::Vector3f color = texture_color;
    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    Eigen::Vector3f result_color = {0, 0, 0};

    for (auto &light: lights) {
        // For each light source in the code, calculate what the *ambient*, *diffuse*, and *specular*
        // components are. Then, accumulate that result on the *result_color* object.

        // blinn-phong模型
        auto n_vec = normal.normalized();
        auto l_vec = (light.position - point).normalized();
        auto v_vec = (eye_pos - point).normalized();
        auto h_vec = (l_vec + v_vec).normalized();
        auto r_dis = std::pow(light.position.x() - point.x(), 2) + std::pow(light.position.y() - point.y(), 2) + std::pow(light.position.z() - point.z(), 2);

        Eigen::Vector3f ambient = ka.cwiseProduct(amb_light_intensity);
        Eigen::Vector3f diffuse = kd.cwiseProduct(light.intensity / r_dis) * std::max(0.0f, n_vec.dot(l_vec));
        Eigen::Vector3f specular = ks.cwiseProduct(light.intensity / r_dis) * std::pow(std::max(0.0f, h_vec.dot(n_vec)), p);
        result_color += ambient + diffuse + specular;
    }

    return result_color * 255.f;
}

Eigen::Vector3f phong_fragment_shader(const fragment_shader_payload &payload) {
    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = payload.color;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20,  20,  20},
                    {500, 500, 500}};
    auto l2 = light{{-20, 20,  0},
                    {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};
    Eigen::Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Eigen::Vector3f color = payload.color;
    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    Eigen::Vector3f result_color = {0, 0, 0};
    for (auto &light: lights) {
        // components are. Then, accumulate that result on the *result_color* object.
        // 公式参考lecture08课件
        auto n_vec = normal.normalized();
        auto l_vec = (light.position - point).normalized();
        auto v_vec = (eye_pos - point).normalized();
        auto h_vec = (l_vec + v_vec).normalized();
        auto r_dis = std::pow(light.position.x() - point.x(), 2) + std::pow(light.position.y() - point.y(), 2) + std::pow(light.position.z() - point.z(), 2);

        Eigen::Vector3f ambient = ka.cwiseProduct(amb_light_intensity);
        Eigen::Vector3f diffuse = kd.cwiseProduct(light.intensity / r_dis) * std::max(0.0f, n_vec.dot(l_vec));
        Eigen::Vector3f specular = ks.cwiseProduct(light.intensity / r_dis) * std::pow(std::max(0.0f, h_vec.dot(n_vec)), p);
        result_color += ambient + diffuse + specular;
    }

    return result_color * 255.f;
}


Eigen::Vector3f displacement_fragment_shader(const fragment_shader_payload &payload) {

    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = payload.color;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20,  20,  20},
                    {500, 500, 500}};
    auto l2 = light{{-20, 20,  0},
                    {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};
    Eigen::Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Eigen::Vector3f color = payload.color;
    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;

    float kh = 0.2, kn = 0.1;

    // Implement displacement mapping here
    // Let n = normal = (x, y, z)
    // Vector t = (x*y/sqrt(x*x+z*z),sqrt(x*x+z*z),z*y/sqrt(x*x+z*z))
    // Vector b = n cross product t
    // Matrix TBN = [t b n]
    // dU = kh * kn * (h(u+1/w,v)-h(u,v))
    // dV = kh * kn * (h(u,v+1/h)-h(u,v))
    // Vector ln = (-dU, -dV, 1)
    // Position p = p + kn * n * h(u,v)
    // Normal n = normalize(TBN * ln)
    Eigen::Vector3f n = normal.normalized();
    Eigen::Vector3f t = {n.x() * n.y() / std::sqrt(n.x() * n.x() + n.z() * n.z()), std::sqrt(n.x() * n.x() + n.z() * n.z()), n.z() * n.y() / std::sqrt(n.x() * n.x() + n.z() * n.z())};
    Eigen::Vector3f b = n.cross(t);
    Eigen::Matrix3f TBN;
    TBN << t, b, n;
    float dU = kh * kn * (payload.texture->getColor(payload.tex_coords.x() + 1.0f / payload.texture->width, payload.tex_coords.y()).norm() - payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y()).norm());
    float dV = kh * kn * (payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y() + 1.0f / payload.texture->height).norm() - payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y()).norm());
    Eigen::Vector3f ln = {-dU, -dV, 1};
    point = payload.view_pos + kn * n * payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y()).norm();
    normal = (TBN * ln).normalized();

    Eigen::Vector3f result_color = {0, 0, 0};

    for (auto &light: lights) {
        // For each light source in the code, calculate what the *ambient*, *diffuse*, and *specular*
        // components are. Then, accumulate that result on the *result_color* object.
        auto n_vec = normal.normalized();
        auto l_vec = (light.position - point).normalized();
        auto v_vec = (eye_pos - point).normalized();
        auto h_vec = (l_vec + v_vec).normalized();
        auto r_dis = std::pow(light.position.x() - point.x(), 2) + std::pow(light.position.y() - point.y(), 2) + std::pow(light.position.z() - point.z(), 2);

        Eigen::Vector3f ambient = ka.cwiseProduct(amb_light_intensity);
        Eigen::Vector3f diffuse = kd.cwiseProduct(light.intensity / r_dis) * std::max(0.0f, n_vec.dot(l_vec));
        Eigen::Vector3f specular = ks.cwiseProduct(light.intensity / r_dis) * std::pow(std::max(0.0f, h_vec.dot(n_vec)), p);
        result_color += ambient + diffuse + specular;
    }

    return result_color * 255.f;
}


Eigen::Vector3f bump_fragment_shader(const fragment_shader_payload &payload) {

    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f kd = payload.color;
    Eigen::Vector3f ks = Eigen::Vector3f(0.7937, 0.7937, 0.7937);

    auto l1 = light{{20,  20,  20},
                    {500, 500, 500}};
    auto l2 = light{{-20, 20,  0},
                    {500, 500, 500}};

    std::vector<light> lights = {l1, l2};
    Eigen::Vector3f amb_light_intensity{10, 10, 10};
    Eigen::Vector3f eye_pos{0, 0, 10};

    float p = 150;

    Eigen::Vector3f color = payload.color;
    Eigen::Vector3f point = payload.view_pos;
    Eigen::Vector3f normal = payload.normal;


    float kh = 0.2, kn = 0.1;

    // Implement bump mapping here
    // Let n = normal = (x, y, z)
    // Vector t = (x*y/sqrt(x*x+z*z),sqrt(x*x+z*z),z*y/sqrt(x*x+z*z))
    // Vector b = n cross product t
    // Matrix TBN = [t b n]
    // dU = kh * kn * (h(u+1/w,v)-h(u,v))
    // dV = kh * kn * (h(u,v+1/h)-h(u,v))
    // Vector ln = (-dU, -dV, 1)
    // Normal n = normalize(TBN * ln)
    Eigen::Vector3f n = normal.normalized();
    Eigen::Vector3f t = {n.x() * n.y() / std::sqrt(n.x() * n.x() + n.z() * n.z()), std::sqrt(n.x() * n.x() + n.z() * n.z()), n.z() * n.y() / std::sqrt(n.x() * n.x() + n.z() * n.z())};
    Eigen::Vector3f b = n.cross(t);
    Eigen::Matrix3f TBN;
    TBN << t, b, n;
    float dU = kh * kn * (payload.texture->getColor(payload.tex_coords.x() + 1.0f / payload.texture->width, payload.tex_coords.y()).norm() - payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y()).norm());
    float dV = kh * kn * (payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y() + 1.0f / payload.texture->height).norm() - payload.texture->getColor(payload.tex_coords.x(), payload.tex_coords.y()).norm());
    Eigen::Vector3f ln = {-dU, -dV, 1};
    normal = (TBN * ln).normalized();

    Eigen::Vector3f result_color = normal;

    return result_color * 255.f;
}

enum SHADER_TYPE {
    SHADER_TYPE_TEXTURE,
    SHADER_TYPE_NORMAL,
    SHADER_TYPE_PHONG,
    SHADER_TYPE_BUMP,
    SHADER_TYPE_DISPLACEMENT
};

struct shader_type_info {
    SHADER_TYPE type;
    std::string texture_path;
    std::function < Eigen::Vector3f(fragment_shader_payload)> active_shader;
    shader_type_info(SHADER_TYPE in_type, std::string in_texture_path, std::function < Eigen::Vector3f(fragment_shader_payload)> in_active_shader)
        : type(in_type), texture_path(std::move(in_texture_path)), active_shader(std::move(in_active_shader)) {}
};

int main(int argc, const char **argv) {
    std::vector<Triangle *> TriangleList;

    float angle = 140.0;
    bool command_line = false;

    std::string working_dir = "../../";
    std::string filename = "output.png";
    objl::Loader Loader;
    std::string obj_path = working_dir + "models/spot/";

    // Load .obj File
    bool loadout = Loader.LoadFile(working_dir + "models/spot/spot_triangulated_good.obj");
    for (auto mesh: Loader.LoadedMeshes) {
        for (int i = 0; i < mesh.Vertices.size(); i += 3) {
            Triangle *t = new Triangle();
            for (int j = 0; j < 3; j++) {
                t->setVertex(j, Vector4f(mesh.Vertices[i + j].Position.X, mesh.Vertices[i + j].Position.Y, mesh.Vertices[i + j].Position.Z, 1.0));
                t->setNormal(j, Vector3f(mesh.Vertices[i + j].Normal.X, mesh.Vertices[i + j].Normal.Y, mesh.Vertices[i + j].Normal.Z));
                t->setTexCoord(j, Vector2f(mesh.Vertices[i + j].TextureCoordinate.X, mesh.Vertices[i + j].TextureCoordinate.Y));
            }
            TriangleList.push_back(t);
        }
    }

    rst::rasterizer r(700, 700);

    // init type
    // change type
    SHADER_TYPE type = SHADER_TYPE_TEXTURE;
    shader_type_info current_shader(SHADER_TYPE_TEXTURE, "spot_texture.png", texture_fragment_shader);

    switch (type) {
        case SHADER_TYPE_TEXTURE:
            current_shader = shader_type_info(SHADER_TYPE_TEXTURE, "spot_texture.png", texture_fragment_shader);
            break;
        case SHADER_TYPE_NORMAL:
            current_shader = shader_type_info(SHADER_TYPE_NORMAL, "hmap.jpg", normal_fragment_shader);
            break;
        case SHADER_TYPE_PHONG:
            current_shader = shader_type_info(SHADER_TYPE_PHONG, "hmap.jpg", phong_fragment_shader);
            break;
        case SHADER_TYPE_BUMP:
            current_shader = shader_type_info(SHADER_TYPE_BUMP, "hmap.jpg", bump_fragment_shader);
            break;
        case SHADER_TYPE_DISPLACEMENT:
            current_shader = shader_type_info(SHADER_TYPE_DISPLACEMENT, "hmap.jpg", displacement_fragment_shader);
            break;
        default:
            break;
    }

    // Change shader type
    auto texture_path = current_shader.texture_path;
    r.set_texture(Texture(obj_path + texture_path));
    std::function < Eigen::Vector3f(fragment_shader_payload) > active_shader = current_shader.active_shader;

    if (argc >= 2) {
        command_line = true;
        filename = std::string(argv[1]);

        if (argc == 3 && std::string(argv[2]) == "texture") {
            std::cout << "Rasterizing using the texture shader\n";
            active_shader = texture_fragment_shader;
            texture_path = "spot_texture.png";
            r.set_texture(Texture(obj_path + texture_path));
        } else if (argc == 3 && std::string(argv[2]) == "normal") {
            std::cout << "Rasterizing using the normal shader\n";
            active_shader = normal_fragment_shader;
        } else if (argc == 3 && std::string(argv[2]) == "phong") {
            std::cout << "Rasterizing using the phong shader\n";
            active_shader = phong_fragment_shader;
        } else if (argc == 3 && std::string(argv[2]) == "bump") {
            std::cout << "Rasterizing using the bump shader\n";
            active_shader = bump_fragment_shader;
        } else if (argc == 3 && std::string(argv[2]) == "displacement") {
            std::cout << "Rasterizing using the bump shader\n";
            active_shader = displacement_fragment_shader;
        }
    }

    Eigen::Vector3f eye_pos = {0, 0, 10};

    r.set_vertex_shader(vertex_shader);
    r.set_fragment_shader(active_shader);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);
        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45.0, 1, 0.1, 50));

        r.draw(TriangleList);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45.0, 1, 0.1, 50));

        //r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);
        r.draw(TriangleList);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imshow("image", image);
        cv::imwrite(filename, image);
        key = cv::waitKey(10);

        if (key == 'a') {
            angle -= 0.1;
        } else if (key == 'd') {
            angle += 0.1;
        }

    }
    return 0;
}
