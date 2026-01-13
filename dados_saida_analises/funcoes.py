
import os
import rasterio
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress 
from rasterio.mask import mask
import geopandas as gpd
import pandas as pd
from scipy.stats import linregress
from pymannkendall import original_test, sens_slope
import numpy as np

def calculo_media_bh_pet(anos, suffix_name, gdf_corte=None):
    medias = []
    anos_validos = []
    estatisticas = {'min_values': [], 'max_values': [], 'std_values': []}
    for ano in anos:
        caminho_tif = f"../input_data/03_et0/ET0_annual/et0_{ano}.tif"


        with rasterio.open(caminho_tif) as src:
            # Recorte por GDF (se fornecido)
            if gdf_corte is not None:
                # Garante que o GDF esteja no mesmo CRS do raster
                gdf_corte = gdf_corte.to_crs(src.crs)
                # Recorta o raster
                dados_recortados, _ = mask(src, gdf_corte.geometry, crop=True, nodata=src.nodata)
                dados = dados_recortados[0]  # Pega a primeira banda
            else:
                dados = src.read(1)
            
            # Processamento dos dados
            nodata = src.nodata
            mask_valido = ~np.isnan(dados)
            if nodata is not None:
                mask_valido = mask_valido & (dados != nodata)
            
            dados_validos = dados[mask_valido]
            
            if len(dados_validos) > 0:
                medias.append(np.mean(dados_validos))
                anos_validos.append(ano)
                estatisticas['min_values'].append(np.min(dados_validos))
                estatisticas['max_values'].append(np.max(dados_validos))
                estatisticas['std_values'].append(np.std(dados_validos))


    medias = np.array(medias)
    estatisticas['min_values'] = np.array(estatisticas['min_values'])
    estatisticas['max_values'] = np.array(estatisticas['max_values'])
    estatisticas['std_values'] = np.array(estatisticas['std_values'])
    
    media_geral = np.mean(medias) if len(medias) > 0 else np.nan
    std_geral = np.mean(estatisticas['std_values']) if len(estatisticas['std_values']) > 0 else np.nan
    estatisticas.update({
        'coef_variacao': (np.std(medias)/media_geral)*100 if len(medias) > 0 else np.nan,
        'anos_processados': anos_validos,
        'n_anos_validos': len(medias)
    })
    
    print(f"Média geral: {media_geral:.2f} mm")
    print(f"Coeficiente de variação: {estatisticas['coef_variacao']:.1f}%")
    print(f"Desvio-padrão geral: {std_geral:.1f}")
    return medias, media_geral, estatisticas


def calculo_media_bh(anos, suffix_name, gdf_corte=None):
    
    medias = []
    anos_validos = []
    estatisticas = {
        'min_values': [],
        'max_values': [],
        'std_values': [],
        'coef_variacao': np.nan,
        'anos_processados': [],
        'n_anos_validos': 0
    }
    
    for ano in anos:
        caminho_tif = os.path.join(
            './final', 
            str(ano), 
            'tfa500_alpha0.0833_beta1_gamma1_flowMFD', 
            f'{suffix_name}_{ano}_tfa500_alpha0.0833_beta1_gamma1_flowMFD.tif'
        )
        
        try:
            with rasterio.open(caminho_tif) as src:
                # Recorte por GDF (se fornecido)
                if gdf_corte is not None:
                    # Garante que o GDF esteja no mesmo CRS do raster
                    gdf_corte = gdf_corte.to_crs(src.crs)
                    # Recorta o raster
                    dados_recortados, _ = mask(src, gdf_corte.geometry, crop=True, nodata=src.nodata)
                    dados = dados_recortados[0]  # Pega a primeira banda
                else:
                    dados = src.read(1)
                
                # Processamento dos dados
                nodata = src.nodata
                mask_valido = ~np.isnan(dados)
                if nodata is not None:
                    mask_valido = mask_valido & (dados != nodata)
                
                dados_validos = dados[mask_valido]
                
                if len(dados_validos) > 0:
                    medias.append(np.mean(dados_validos))
                    anos_validos.append(ano)
                    estatisticas['min_values'].append(np.min(dados_validos))
                    estatisticas['max_values'].append(np.max(dados_validos))
                    estatisticas['std_values'].append(np.std(dados_validos))
                    
        except Exception as e:
            print(f"Erro no ano {ano}: {str(e)}")
            continue
    
    # Conversão para arrays numpy
    medias = np.array(medias)
    estatisticas['min_values'] = np.array(estatisticas['min_values'])
    estatisticas['max_values'] = np.array(estatisticas['max_values'])
    estatisticas['std_values'] = np.array(estatisticas['std_values'])
    
    # Cálculo das estatísticas gerais
    media_geral = np.mean(medias) if len(medias) > 0 else np.nan
    std_geral = np.mean(estatisticas['std_values']) if len(estatisticas['std_values']) > 0 else np.nan
    
    estatisticas.update({
        'coef_variacao': (np.std(medias)/media_geral)*100 if len(medias) > 0 and media_geral != 0 else np.nan,
        'anos_processados': anos_validos,
        'n_anos_validos': len(medias)
    })
    
    print(f"\nResultados para {suffix_name}:")
    print(f"- Média geral: {media_geral:.2f} mm")
    print(f"- Coeficiente de variação: {estatisticas['coef_variacao']:.1f}%")
    print(f"- Desvio-padrão geral: {std_geral:.1f}")
    
    return medias, media_geral, estatisticas
def analise_tendencia(anos, medias, local="BH"):
    anos = np.array(anos)
    medias = np.array(medias)
    # --- 1. Regressão Linear ---
    slope, intercept, r_value, p_value, _ = linregress(anos, medias)
    r2 = r_value**2
    tendencia = intercept + slope * anos  
    # --- 2. Testes Não-Paramétricos ---
    try:
        mk_result = original_test(medias)
        sen_slope, sen_intercept = sens_slope(medias)
        mk_trend = mk_result.trend
        mk_p = mk_result.p
    except ImportError:
        mk_trend = "N/A (instale pymannkendall)"
        mk_p = np.nan
        sen_slope = np.nan   
    # --- 3. Análise por Subperíodos ---
    def analisar_subperiodo(mask, nome):
        if sum(mask) < 5:  # Mínimo de 5 anos para análise
            return {f'slope_{nome}': np.nan, f'p_{nome}': np.nan}
        slope, _, _, p, _ = linregress(anos[mask], medias[mask])
        return {f'slope_{nome}': slope, f'p_{nome}': p}
    # Análise antes/depois de 2000
    pre2000 = anos < 2000
    post2000 = anos >= 2000
    subperiodos = {
        **analisar_subperiodo(pre2000, 'pre2000'),
        **analisar_subperiodo(post2000, 'post2000')
    }
    # --- 4. Compilar Resultados ---
    resultados = {
        'local': local,
        'n_anos': len(medias),
        'periodo': f"{min(anos)}-{max(anos)}",
        'media_geral': np.mean(medias),
        'slope': slope,
        'intercept': intercept,
        'p_value': p_value,
        'r2': r2,
        'mk_trend': mk_trend,
        'mk_p': mk_p,
        'sen_slope': sen_slope,
        **subperiodos,
        'dados': {
            'anos': anos,
            'medias': medias,
            'tendencia': tendencia
        }
    }   
    # --- 5. Relatório Automático ---
    print(f"Análise de Tendência para {local} ({resultados['periodo']})")
    print("----------------------------------------")
    print("Regressão Linear:")
    print(f" • Inclinação: {slope:.2f} mm/ano ({'↑ aumento' if slope > 0 else '↓ redução'})")
    print(f" • p-valor: {p_value:.3f} {'(significativo)' if p_value < 0.05 else '(não significativo)'}")
    print(f" • R²: {r2:.2f} (explica {r2*100:.0f}% da variabilidade)")
    
    if 'sen_slope' in resultados:
        print("Teste de Mann-Kendall:")
        print(f" • Tendência: {mk_trend} (p={mk_p:.3f})")
        print(f" • Sen's Slope: {sen_slope:.2f} mm/ano")
    
    print("Análise por Subperíodos:")
    print(f" • 1985-1999: {resultados['slope_pre2000']:.2f} mm/ano (p={resultados['p_pre2000']:.3f})")
    print(f" • 2000-2023: {resultados['slope_post2000']:.2f} mm/ano (p={resultados['p_post2000']:.3f})")
    
    return resultados


def plot_medias_bh_simples(resultados, variavel, suffix_name):
    plt.figure(figsize=(14, 7))
    # 1. Dados anuais
    plt.plot(resultados['dados']['anos'], resultados['dados']['medias'], 
            'o-', color='#3498db', markersize=5, linewidth=1.5, 
            label='Dados anuais', alpha=0.8)

    # 4. Configurações do gráfico
    plt.title(f'{variavel} na Bacia do Paranoá ({resultados["periodo"]})', 
            pad=20, fontsize=14)
    plt.xlabel('Ano', fontsize=12)
    plt.ylabel(f'{variavel} média (mm)', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.3)
        
    # 6. Salvar em alta resolução (escolha UM dos formatos)
    plt.savefig(f'./graficos/{suffix_name.split("./intermediate_outputs")[-1]}_paranoa_simples.png', dpi=300, bbox_inches='tight')

    plt.tight_layout()
    plt.show()

def plot_medias_bh(resultados, variavel, suffix_name):
    plt.figure(figsize=(14, 7))
    # 1. Dados anuais
    plt.plot(resultados['dados']['anos'], resultados['dados']['medias'], 
            'o-', color='#3498db', markersize=5, linewidth=1.5, 
            label='Média anual', alpha=0.8)

    # 2. Linha de tendência (inclinação de -4.3 mm/ano)
    plt.plot(resultados['dados']['anos'], resultados['dados']['tendencia'], 
            '--', color='#e74c3c', linewidth=1.8,
            label=f'Tendência ({resultados["slope"]:.1f} mm/ano)')

    # 3. Média geral
    plt.axhline(y=resultados['media_geral'], color='#2ecc71', 
                linestyle=':', linewidth=1.8,
                label=f'Média geral ({resultados["media_geral"]:.1f} mm)')

    # 4. Configurações do gráfico
    plt.title(f'{variavel} na Bacia do Paranoá ({resultados["periodo"]})', 
            pad=20, fontsize=14)
    plt.xlabel('Ano', fontsize=12)
    plt.ylabel(f'{variavel} (mm)', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.3)

    # 5. Legenda em quadro no canto superior direito
    leg = plt.legend(frameon=True, fontsize=12, 
                    loc= "upper right", #'lower right',
                    bbox_to_anchor=(0.98, 0.98), #(0.98, 0.02), #
                    edgecolor='gray',              
                    shadow=False,               
                    framealpha=0.9)         
    # 6. Salvar em alta resolução (escolha UM dos formatos)
    plt.savefig(f'./graficos/{suffix_name.split("./intermediate_outputs")[-1]}_paranoa_completo.png', dpi=1200, bbox_inches='tight')

    plt.tight_layout()
    plt.show()
def plot_medias_bh_lower(resultados, variavel, suffix_name):
    plt.figure(figsize=(14, 7))
    # 1. Dados anuais
    plt.plot(resultados['dados']['anos'], resultados['dados']['medias'], 
            'o-', color='#3498db', markersize=5, linewidth=1.5, 
            label='Média anual', alpha=0.8)

    # 2. Linha de tendência (inclinação de -4.3 mm/ano)
    plt.plot(resultados['dados']['anos'], resultados['dados']['tendencia'], 
            '--', color='#e74c3c', linewidth=1.8,
            label=f'Tendência ({resultados["slope"]:.1f} mm/ano)')

    # 3. Média geral
    plt.axhline(y=resultados['media_geral'], color='#2ecc71', 
                linestyle=':', linewidth=1.8,
                label=f'Média geral ({resultados["media_geral"]:.1f} mm)')

    # 4. Configurações do gráfico
    plt.title(f'{variavel} na Bacia do Paranoá ({resultados["periodo"]})', 
            pad=20, fontsize=14)
    plt.xlabel('Ano', fontsize=12)
    plt.ylabel(f'{variavel} (mm)', fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.3)

    # 5. Legenda em quadro no canto superior direito
    leg = plt.legend(frameon=True, fontsize=12, 
                    loc= 'lower right',
                    bbox_to_anchor=(0.98, 0.02), #
                    edgecolor='gray',              
                    shadow=False,               
                    framealpha=0.9)         
    # 6. Salvar em alta resolução (escolha UM dos formatos)
    plt.savefig(f'./graficos/{suffix_name.split("./intermediate_outputs")[-1]}_paranoa_completo.png', dpi=1200, bbox_inches='tight')

    plt.tight_layout()
    plt.show()


def calcular_medias_zonas_paranoa(gdf, anos, suffix_name, zona_column='uh_label', verbose=True):

    resultados = {
        'zonas': gdf[zona_column].unique().tolist(),
        'dados': {},
        'estatisticas': {}
    }
    
    # Listas para armazenar dados das tabelas
    dados_estatisticas = []
    dados_tendencias = []

    for _, zona in gdf.iterrows():
        zona_id = zona[zona_column]
        geometria = zona['geometry']
        
        medias = []
        estatisticas = {
            'min_values': [],
            'max_values': [],
            'std_values': [],
            'anos_faltantes': []
        }

        for ano in anos:
            caminho_tif = os.path.join('./final', str(ano), 
                                    'tfa500_alpha0.0833_beta1_gamma1_flowMFD', 
                                    f'{suffix_name}_{ano}_tfa500_alpha0.0833_beta1_gamma1_flowMFD.tif')
            
            if not os.path.exists(caminho_tif):
                estatisticas['anos_faltantes'].append(ano)
                continue

            try:
                with rasterio.open(caminho_tif) as src:
                    dados_recortados, _ = mask(src, [geometria], crop=True, nodata=src.nodata)
                    dados = dados_recortados[0]
                    
                    dados_validos = dados[~np.isnan(dados)]
                    if src.nodata is not None:
                        dados_validos = dados_validos[dados_validos != src.nodata]
                    
                    if len(dados_validos) == 0:
                        estatisticas['anos_faltantes'].append(ano)
                        continue
                    
                    medias.append(np.mean(dados_validos))
                    estatisticas['min_values'].append(np.min(dados_validos))
                    estatisticas['max_values'].append(np.max(dados_validos))
                    estatisticas['std_values'].append(np.std(dados_validos))
            
            except Exception as e:
                estatisticas['anos_faltantes'].append(ano)
                continue

        # Processamento final
        medias = np.array(medias)
        estatisticas['min_values'] = np.array(estatisticas['min_values'])
        estatisticas['max_values'] = np.array(estatisticas['max_values'])
        estatisticas['std_values'] = np.array(estatisticas['std_values'])
        
        media_geral = np.nanmean(medias)
        std_geral = np.nanmean(estatisticas['std_values'])
        cv = (np.nanstd(medias)/media_geral)*100 if media_geral > 0 else np.nan
        
        estatisticas.update({
            'media_geral': media_geral,
            'std_geral': std_geral,
            'coef_variacao': cv,
            'n_anos_validos': len(medias)
        })

        # Armazena resultados brutos
        resultados['dados'][zona_id] = {
            'medias': medias,
            'anos_validos': [ano for ano in anos if ano not in estatisticas['anos_faltantes']]
        }
        resultados['estatisticas'][zona_id] = estatisticas

        # Prepara dados para tabela estatística
        dados_estatisticas.append({
            'Zona': zona_id,
            'Média Geral (mm)': f"{media_geral:.2f}",
            'CV (%)': f"{cv:.1f}",
            'Variação Anual (±mm)': f"{std_geral:.1f}",
            'Anos Processados': f"{len(medias)}/{len(anos)}"
        })

        # Calcula tendências para a tabela
        slope, _, r_value, p_value, _ = linregress(range(len(medias)), medias)
        try:
            mk_test = original_test(medias)
            sen_slope, _ = sens_slope(medias)
            mk_trend = mk_test.trend
            mk_p = mk_test.p
        except:
            sen_slope = np.nan
            mk_trend = "N/A"
            mk_p = np.nan

        dados_tendencias.append({
            'Zona': zona_id,
            'Tendência (mm/ano)': f"{slope:.2f}",
            'p-valor': f"{p_value:.3f}",
            'R²': f"{r_value**2:.2f}",
            "Sen's Slope": f"{sen_slope:.2f}",
            'Mann-Kendall': f"{mk_trend} (p={mk_p:.3f})"
        })

    # Cria DataFrames
    resultados['tabela_estatisticas'] = pd.DataFrame(dados_estatisticas)
    resultados['tabela_tendencias'] = pd.DataFrame(dados_tendencias)

    if verbose:
        print("\n=== Estatísticas Descritivas ===")
        print(resultados['tabela_estatisticas'].to_markdown(index=False))
        
        print("\n=== Análise de Tendências ===")
        print(resultados['tabela_tendencias'].to_markdown(index=False))

    return resultados


def calcular_medias_zonas_paranoa_pet(gdf, anos, zona_column='uh_label', verbose=True):

    resultados = {
        'zonas': gdf[zona_column].unique().tolist(),
        'dados': {},
        'estatisticas': {}
    }
    
    # Listas para armazenar dados das tabelas
    dados_estatisticas = []
    dados_tendencias = []

    for _, zona in gdf.iterrows():
        zona_id = zona[zona_column]
        geometria = zona['geometry']
        
        medias = []
        estatisticas = {
            'min_values': [],
            'max_values': [],
            'std_values': [],
            'anos_faltantes': []
        }

        for ano in anos:
            caminho_tif = f"../input_data/03_et0/ET0_annual/et0_{ano}.tif"
            
            if not os.path.exists(caminho_tif):
                estatisticas['anos_faltantes'].append(ano)
                continue

            try:
                with rasterio.open(caminho_tif) as src:
                    dados_recortados, _ = mask(src, [geometria], crop=True, nodata=src.nodata)
                    dados = dados_recortados[0]
                    
                    dados_validos = dados[~np.isnan(dados)]
                    if src.nodata is not None:
                        dados_validos = dados_validos[dados_validos != src.nodata]
                    
                    if len(dados_validos) == 0:
                        estatisticas['anos_faltantes'].append(ano)
                        continue
                    
                    medias.append(np.mean(dados_validos))
                    estatisticas['min_values'].append(np.min(dados_validos))
                    estatisticas['max_values'].append(np.max(dados_validos))
                    estatisticas['std_values'].append(np.std(dados_validos))
            
            except Exception as e:
                estatisticas['anos_faltantes'].append(ano)
                continue

        # Processamento final
        medias = np.array(medias)
        estatisticas['min_values'] = np.array(estatisticas['min_values'])
        estatisticas['max_values'] = np.array(estatisticas['max_values'])
        estatisticas['std_values'] = np.array(estatisticas['std_values'])
        
        media_geral = np.nanmean(medias)
        std_geral = np.nanmean(estatisticas['std_values'])
        cv = (np.nanstd(medias)/media_geral)*100 if media_geral > 0 else np.nan
        
        estatisticas.update({
            'media_geral': media_geral,
            'std_geral': std_geral,
            'coef_variacao': cv,
            'n_anos_validos': len(medias)
        })

        # Armazena resultados brutos
        resultados['dados'][zona_id] = {
            'medias': medias,
            'anos_validos': [ano for ano in anos if ano not in estatisticas['anos_faltantes']]
        }
        resultados['estatisticas'][zona_id] = estatisticas

        # Prepara dados para tabela estatística
        dados_estatisticas.append({
            'Zona': zona_id,
            'Média Geral (mm)': f"{media_geral:.2f}",
            'CV (%)': f"{cv:.1f}",
            'Variação Anual (±mm)': f"{std_geral:.1f}",
            'Anos Processados': f"{len(medias)}/{len(anos)}"
        })

        # Calcula tendências para a tabela
        slope, _, r_value, p_value, _ = linregress(range(len(medias)), medias)
        try:
            mk_test = original_test(medias)
            sen_slope, _ = sens_slope(medias)
            mk_trend = mk_test.trend
            mk_p = mk_test.p
        except:
            sen_slope = np.nan
            mk_trend = "N/A"
            mk_p = np.nan

        dados_tendencias.append({
            'Zona': zona_id,
            'Tendência (mm/ano)': f"{slope:.2f}",
            'p-valor': f"{p_value:.3f}",
            'R²': f"{r_value**2:.2f}",
            "Sen's Slope": f"{sen_slope:.2f}",
            'Mann-Kendall': f"{mk_trend} (p={mk_p:.3f})"
        })

    # Cria DataFrames
    resultados['tabela_estatisticas'] = pd.DataFrame(dados_estatisticas)
    resultados['tabela_tendencias'] = pd.DataFrame(dados_tendencias)

    if verbose:
        print("\n=== Estatísticas Descritivas ===")
        print(resultados['tabela_estatisticas'].to_markdown(index=False))
        
        print("\n=== Análise de Tendências ===")
        print(resultados['tabela_tendencias'].to_markdown(index=False))

    return resultados

       
def calcular_medias_zonas_paranoa_temp(gdf, anos, suffix_name, zona_column='uh_label', verbose=False):
   
    resultados = {
        'zonas': gdf[zona_column].unique().tolist(),
        'dados': {},
        'estatisticas': {}
    }

    for _, zona in gdf.iterrows():
        zona_id = zona[zona_column]
        geometria = zona['geometry']
        
        medias = []
        estatisticas = {
            'min_values': [],
            'max_values': [],
            'std_values': [],
            'anos_faltantes': []
        }

        for ano in anos:
            caminho_tif = os.path.join('./final', str(ano), 
                                    'tfa500_alpha0.0833_beta1_gamma1_flowMFD', 
                                    f'{suffix_name}_{ano}_tfa500_alpha0.0833_beta1_gamma1_flowMFD.tif')
            
            if not os.path.exists(caminho_tif):
                estatisticas['anos_faltantes'].append(ano)
                continue

            try:
                with rasterio.open(caminho_tif) as src:
                    # Recorta para a zona
                    dados_recortados, _ = mask(src, [geometria], crop=True, nodata=src.nodata)
                    dados = dados_recortados[0]
                    
                    # Filtra valores válidos
                    dados_validos = dados[~np.isnan(dados)]
                    if src.nodata is not None:
                        dados_validos = dados_validos[dados_validos != src.nodata]
                    
                    if len(dados_validos) == 0:
                        estatisticas['anos_faltantes'].append(ano)
                        continue
                    
                    # Calcula estatísticas
                    medias.append(np.mean(dados_validos))
                    estatisticas['min_values'].append(np.min(dados_validos))
                    estatisticas['max_values'].append(np.max(dados_validos))
                    estatisticas['std_values'].append(np.std(dados_validos))
            
            except Exception as e:
                if verbose:
                    print(f"\nErro na zona {zona_id}, ano {ano}: {str(e)}")
                estatisticas['anos_faltantes'].append(ano)
                continue

        # Converte para arrays numpy
        medias = np.array(medias)
        estatisticas['min_values'] = np.array(estatisticas['min_values'])
        estatisticas['max_values'] = np.array(estatisticas['max_values'])
        estatisticas['std_values'] = np.array(estatisticas['std_values'])
        
        # Estatísticas consolidadas
        media_geral = np.nanmean(medias) if len(medias) > 0 else np.nan
        std_geral = np.nanmean(estatisticas['std_values']) if len(estatisticas['std_values']) > 0 else np.nan
        
        estatisticas.update({
            'media_geral': media_geral,
            'std_geral': std_geral,
            'coef_variacao': (np.nanstd(medias)/media_geral)*100 if media_geral > 0 else np.nan,
            'n_anos_validos': len(medias)
        })

        # Armazena resultados
        resultados['dados'][zona_id] = {
            'medias': medias,
            'anos_validos': [ano for ano in anos if ano not in estatisticas['anos_faltantes']]
        }
        
        resultados['estatisticas'][zona_id] = estatisticas


    return resultados

def calcular_medias_zonas_paranoa_temp_pet(gdf, anos, zona_column='uh_label', verbose=False):
   
    resultados = {
        'zonas': gdf[zona_column].unique().tolist(),
        'dados': {},
        'estatisticas': {}
    }

    for _, zona in gdf.iterrows():
        zona_id = zona[zona_column]
        geometria = zona['geometry']
        
        medias = []
        estatisticas = {
            'min_values': [],
            'max_values': [],
            'std_values': [],
            'anos_faltantes': []
        }

        for ano in anos:
            caminho_tif = f"../input_data/03_et0/ET0_annual/et0_{ano}.tif"
            if not os.path.exists(caminho_tif):
                estatisticas['anos_faltantes'].append(ano)
                continue

            try:
                with rasterio.open(caminho_tif) as src:
                    # Recorta para a zona
                    dados_recortados, _ = mask(src, [geometria], crop=True, nodata=src.nodata)
                    dados = dados_recortados[0]
                    
                    # Filtra valores válidos
                    dados_validos = dados[~np.isnan(dados)]
                    if src.nodata is not None:
                        dados_validos = dados_validos[dados_validos != src.nodata]
                    
                    if len(dados_validos) == 0:
                        estatisticas['anos_faltantes'].append(ano)
                        continue
                    
                    # Calcula estatísticas
                    medias.append(np.mean(dados_validos))
                    estatisticas['min_values'].append(np.min(dados_validos))
                    estatisticas['max_values'].append(np.max(dados_validos))
                    estatisticas['std_values'].append(np.std(dados_validos))
            
            except Exception as e:
                if verbose:
                    print(f"\nErro na zona {zona_id}, ano {ano}: {str(e)}")
                estatisticas['anos_faltantes'].append(ano)
                continue

        # Converte para arrays numpy
        medias = np.array(medias)
        estatisticas['min_values'] = np.array(estatisticas['min_values'])
        estatisticas['max_values'] = np.array(estatisticas['max_values'])
        estatisticas['std_values'] = np.array(estatisticas['std_values'])
        
        # Estatísticas consolidadas
        media_geral = np.nanmean(medias) if len(medias) > 0 else np.nan
        std_geral = np.nanmean(estatisticas['std_values']) if len(estatisticas['std_values']) > 0 else np.nan
        
        estatisticas.update({
            'media_geral': media_geral,
            'std_geral': std_geral,
            'coef_variacao': (np.nanstd(medias)/media_geral)*100 if media_geral > 0 else np.nan,
            'n_anos_validos': len(medias)
        })

        # Armazena resultados
        resultados['dados'][zona_id] = {
            'medias': medias,
            'anos_validos': [ano for ano in anos if ano not in estatisticas['anos_faltantes']]
        }
        
        resultados['estatisticas'][zona_id] = estatisticas


    return resultados
    
    
def analise_tendencia_por_zona(resultados, zona_column='uh_label'):

    tendencias = {}
    
    for zona in resultados['zonas']:
        anos_validos = np.array(resultados['dados'][zona]['anos_validos'])
        medias = resultados['dados'][zona]['medias']
        
        # Regressão Linear
        slope, intercept, r_value, p_value, _ = linregress(range(len(anos_validos)), medias)
        r2 = r_value**2
        
        # Teste de Mann-Kendall
        try:
            mk_result = original_test(medias)
            sen_slope_val, _ = sens_slope(medias)
            mk_trend = mk_result.trend
            mk_p = mk_result.p
        except:
            mk_trend = "N/A"
            mk_p = np.nan
            sen_slope_val = np.nan
        
        # Armazenar resultados
        tendencias[zona] = {
            'slope': slope,
            'intercept': intercept,
            'p_value': p_value,
            'r2': r2,
            'mk_trend': mk_trend,
            'mk_p': mk_p,
            'sen_slope': sen_slope_val,
            'anos_validos': anos_validos,
            'medias': medias,
            'tendencia': intercept + slope * np.arange(len(anos_validos))
        }
        

    return tendencias

def plot_tendencias_zonas(tendencias, variavel, suffix_name, zonas=None, figsize=(34, 20)):

    # Configuração inicial
    zonas = zonas or list(tendencias.keys())
    n_zonas = len(zonas)
    n_cols = 2
    n_rows = (n_zonas + 1) // n_cols
    
    # Cria figura e eixos
    fig, axs = plt.subplots(n_rows, n_cols, figsize=figsize)
    axs = axs.flatten() if isinstance(axs, np.ndarray) else [axs]
    
    # Cores e configurações visuais
    colors = {
        'data': '#3498db',
        'trend': '#e74c3c',
        'mean': '#2ecc71'
    }
    line_styles = {
        'data': 'o-',
        'trend': '--',
        'mean': ':'
    }
    
    # Plota dados para cada zona
    for idx, zona in enumerate(zonas):
        ax = axs[idx]
        dados = tendencias[zona]
        
        # Plot dos elementos
        ax.plot(dados['anos_validos'], dados['medias'], 
                line_styles['data'], color=colors['data'], 
                markersize=5, label='Média anual')
        
        ax.plot(dados['anos_validos'], dados['tendencia'], 
                line_styles['trend'], color=colors['trend'], 
                linewidth=1.8, label=f'Tendência ({dados["slope"]:.1f} mm/ano)')
        
        ax.axhline(y=np.mean(dados['medias']), 
                   color=colors['mean'], linestyle=line_styles['mean'],
                   linewidth=1.8, label=f'Média geral ({np.mean(dados["medias"]):.1f} mm)')
        
        # Configurações do gráfico
        ax.set_title(f'UH {zona}', fontsize=28, pad=10)
        ax.set_xlabel('Ano', fontsize=24)
        ax.set_ylabel(f'{variavel} (mm)', fontsize=24)
        ax.tick_params(axis='both', which='major', labelsize=22) 
        ax.grid(axis='y', linestyle='--', alpha=0.3)
        ax.legend(fontsize=22)
               #   loc= "lower right", bbox_to_anchor=(0.98, 0.02))
    
    # Desativa eixos não utilizados
    for j in range(len(zonas), len(axs)):
        axs[j].axis('off')
    
    # Ajustes finais e salvamento
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.33, wspace=0.18)
    
    output_path = f'./graficos/{suffix_name.split("./intermediate_outputs")[-1]}_paranoa_uh_completo.png'
    plt.savefig(output_path, dpi=1200, bbox_inches='tight')
    plt.show()


def plot_tendencias_zonas_simples(tendencias, variavel, suffix_name, zonas=None, figsize=(18, 12)):

    if zonas is None:
        zonas = list(tendencias.keys())
    
    n_zonas = len(zonas)
    n_cols = 2
    n_rows = (n_zonas + 1) // n_cols
    
    fig, axs = plt.subplots(n_rows, n_cols, figsize=figsize)
    axs = axs.flatten() if n_zonas > 1 else [axs]
    
    for idx, zona in enumerate(zonas):
        
        ax = axs[idx]
        dados = tendencias[zona]
        ax.plot(dados['anos_validos'], dados['medias'], 'o-', 
                color='#3498db', markersize=5, label='Média anual')
        ax.set_title(f'Zona {zona}', fontsize=13, pad=12)
        ax.set_xlabel('Ano', fontsize=13)
        ax.set_ylabel(f'{variavel} (mm)', fontsize=13)
        ax.grid(axis='y', linestyle='--', alpha=0.3)       
        ax.legend(fontsize=12)

    if n_zonas % n_cols != 0:
        for idx in range(n_zonas, n_rows * n_cols):
            axs[idx].set_visible(False)  # Desativa os eixos não usados

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3, wspace=0.15)
    plt.savefig(f'./graficos/{suffix_name.split("./intermediate_outputs")[-1]}_paranoa_uh_simples.png', dpi=1500, bbox_inches='tight')
    plt.show()

