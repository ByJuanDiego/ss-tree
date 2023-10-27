import torch
from timm.models import create_model
from torchvision.datasets import ImageFolder
from torch.utils.data import DataLoader
from tqdm import tqdm
import h5py

INPUT_SIZE = 224
BATCH_SIZE = 64
WEIGHT_PATH = "/media/victor/Documentos/Artificio/Demo/artifacts/checkpoint_epoch5.pt"

def torchvision_preprocess_input(image_size, **kwargs):
    from torchvision import transforms
    return transforms.Compose([
        transforms.Resize((INPUT_SIZE, INPUT_SIZE)),
        transforms.CenterCrop((image_size,image_size)),
        torchvision_preprocess(**kwargs),
    ])

def torchvision_preprocess(normalize_mean=(0.485, 0.456, 0.406), normalize_std=(0.229, 0.224, 0.225)):
    from torchvision import transforms
    return transforms.Compose([
        transforms.ToTensor(),
        transforms.Normalize(mean=normalize_mean, std=normalize_std)
    ])

net = create_model('crossvit_18_dagger_408',pretrained = False)
net.cuda()
checkpoint = torch.load(WEIGHT_PATH,map_location=torch.device('cuda:0'))
net.load_state_dict(checkpoint['state_dict'],strict = True)
net.eval()

def get_image_loader(root_dir, transform, batch_size):
    dataset = ImageFolder(root=root_dir, transform=transform)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=4)  # pin_memory=True if using CUDA
    return loader

root_dir = "img"
transform = torchvision_preprocess_input(INPUT_SIZE)
image_loader = get_image_loader(root_dir, transform, BATCH_SIZE)

with h5py.File("embedding.hdf5", "w") as f:
    for batch_idx, (inputs, labels) in tqdm(enumerate(image_loader)):
        inputs = inputs.cuda()
        
        # Ejecutar en modo evaluacion
        with torch.no_grad():
            outputs = net(inputs)
            outputs = outputs.cpu()
        # Suponiendo que outputs es un tensor con los embeddings
        for i, output in enumerate(outputs):
            # Aquí supondré que dataset.samples[i][0] te da la ruta de la imagen
            img_path = image_loader.dataset.samples[batch_idx * BATCH_SIZE + i][0]
            # Creas un grupo por cada imagen
            grp = f.create_group(img_path)
            # Almacenas el embedding y la ruta en ese grupo
            grp.create_dataset("embedding", data=output.numpy())
            grp.attrs["path"] = img_path

            print(img_path)
