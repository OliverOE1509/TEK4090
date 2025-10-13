from controller import Robot

TIME_STEP = 32
robot = Robot()

# Get motors
front_left_motor  = robot.getDevice('front left propeller')
front_right_motor = robot.getDevice('front right propeller')
rear_left_motor   = robot.getDevice('rear left propeller')
rear_right_motor  = robot.getDevice('rear right propeller')

motors = [front_left_motor, front_right_motor, rear_left_motor, rear_right_motor]

# Enable all motors in velocity mode
for m in motors:
    m.setPosition(float('inf'))
    m.setVelocity(0.0)

t = 0.0
while robot.step(TIME_STEP) != -1:
    t += TIME_STEP / 1000.0  # seconds

    # --- Takeoff → hover → land sequence ---
    if t < 2.0:       # takeoff
        vel = 800.0
    elif t < 6.0:     # hover
        vel = 720.0
    elif t < 8.0:     # landing
        vel = 600.0
    else:             # stop
        vel = 0.0

    for m in motors:
        m.setVelocity(vel)
