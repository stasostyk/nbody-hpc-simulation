#pragma once

class RenderBufferObject
{
public:
    virtual void initialize(int n) = 0;
    virtual void draw(int n) = 0;
};
