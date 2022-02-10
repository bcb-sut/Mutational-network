import matplotlib.pyplot as plt


def show_multiple_imgs(imgs, rows, cols):
    fig = plt.figure()
    cnt = 1
    for img in imgs:
        fig.add_subplot(rows, cols, cnt)
        plt.imshow(img)
        plt.axis('off')
        cnt += 1
    plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)


class VisUtils:
    def __init__(self, save_path):
        self.save_path = save_path

    def save_plotly(self, fig, name, frmt):
        fig.write_image(f'{self.save_path}/{frmt}/{name}.{frmt}')

    def save_and_show_plotly(self, fig, name):
        self.save_plotly(fig, name, 'png')
        self.save_plotly(fig, name, 'svg')
        # fig.show()

    def save_plt_with_format(self, name, frmt):
        plt.savefig(f'{self.save_path}/{frmt}/{name}.{frmt}')

    def save_and_show_plt(self, name):
        self.save_plt_with_format(name, 'svg')
        self.save_plt_with_format(name, 'png')

